"""
ProcessRunner - 外部程序进程管理器

提供实时输出流处理、资源监控、超时处理等功能
"""

import os
import shutil
import signal
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from types import FrameType
from typing import IO, Any, Callable, Dict, List, Optional

from ..core.exceptions import ExecutionError, ExecutionTimeoutError
from .logging import get_logger


@dataclass
class ProcessResult:
    """进程执行结果"""

    returncode: int
    stdout: str
    stderr: str
    execution_time: float
    peak_memory_mb: Optional[float] = None


class _OutputBudget:
    """跨 stdout/stderr 读线程共享的输出体积预算（审阅项 C-34）。

    ``run()`` 用两个线程分别抽干 stdout 与 stderr，而 ``max_output_size``
    是防止无界输出把父进程撑爆的护栏。原实现让两个线程对同一个 nonlocal
    整数做裸 ``total += n``：该语句在字节码层面是 LOAD → BINARY_OP → STORE
    三步，语言规范不保证其间不被另一个线程插入，被覆盖的那次增量恰好就是
    "是否已经越线"的依据。

    需要说清楚风险量级：在 GIL 打开的 CPython 上，解释器的特化层会把
    LOAD/BINARY_OP/STORE 融合成不带 eval-breaker 检查的 macro op，因此实测
    很难重现丢更新（本项目 3.14 环境即如此，见 tests/unit/infrastructure/
    test_process_runner.py 的说明）；但 3.13/3.14 起的自由线程构建（PEP 703）
    不再提供这层意外保护。护栏本身是防 OOM 的，不该把正确性押在解释器的
    指令融合细节上，所以这里把"累加 + 越线判定"合成一把锁下的单个原子操作：

    - 每一次累加都对另一个线程可见，计数不丢；
    - 越线只被判定一次，由**先**越线的那个线程拿到 ``True`` 去 terminate；
    - ``error`` 只会被写入一次（原实现里两线程各写一次、后写者胜出）。

    锁的临界区只有一次整数加法与两次比较，纳秒级，不构成吞吐瓶颈。
    """

    __slots__ = ("_max_size", "_lock", "_total", "_exceeded", "_error")

    def __init__(self, max_size: int) -> None:
        self._max_size = max_size
        self._lock = threading.Lock()
        self._total = 0
        self._exceeded = threading.Event()
        self._error: Optional[ExecutionError] = None

    def record(self, size: int) -> bool:
        """累加 ``size`` 字节；返回 True 表示本次累加使总量越过上限。"""
        with self._lock:
            self._total += size
            if self._exceeded.is_set() or self._total <= self._max_size:
                return False
            self._error = ExecutionError(
                f"Output size exceeded maximum limit ({self._max_size} bytes)"
            )
            self._exceeded.set()
            return True

    @property
    def exceeded(self) -> bool:
        """是否已越过上限（任一线程触发即为 True，可安全跨线程读取）。"""
        return self._exceeded.is_set()

    @property
    def total(self) -> int:
        """当前累计字节数（快照值，仅供日志/诊断使用）。"""
        with self._lock:
            return self._total

    @property
    def error(self) -> Optional[ExecutionError]:
        with self._lock:
            return self._error


class ProcessRunner:
    """
    外部程序运行管理器

    特性：
    - 实时输出流处理
    - 内存峰值监控（可选 psutil）
    - 超时处理
    - 日志记录
    - 输出大小限制（stdout/stderr 共享同一把锁保护的预算计数）
    - 可选的返回码强制检查（``run(check=True)``）
    - 信号处理（SIGTERM/SIGINT）
    - macOS/Linux 动态库路径自动处理
    - 优化的I/O缓冲（支持大文件高效处理）
    """

    MAX_OUTPUT_SIZE = 100 * 1024 * 1024
    # 使用 256KB 缓冲区优化大文件I/O性能
    BUFFER_SIZE = 256 * 1024

    def __init__(
        self,
        timeout: Optional[int] = None,
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None,
        monitor_memory: bool = True,
        max_output_size: Optional[int] = None,
        buffer_size: Optional[int] = None,
    ) -> None:
        self.timeout = timeout
        self.cwd = cwd
        self.env = self._prepare_environment(env)
        self.monitor_memory = monitor_memory
        self.max_output_size = max_output_size or self.MAX_OUTPUT_SIZE
        self.buffer_size = buffer_size or self.BUFFER_SIZE
        self.logger = get_logger()
        self._peak_memory = 0.0
        self._memory_monitor_thread: Optional[threading.Thread] = None
        self._stop_monitoring: Optional[threading.Event] = None
        self._process: Optional[subprocess.Popen] = None
        self._original_sigterm_handler: Optional[Any] = None
        self._original_sigint_handler: Optional[Any] = None

    def _prepare_environment(
        self, env: Optional[Dict[str, str]]
    ) -> Optional[Dict[str, str]]:
        """
        准备执行环境，自动处理 macOS 和 Linux 之间的动态链接库路径差异

        macOS 使用 DYLD_LIBRARY_PATH，Linux 使用 LD_LIBRARY_PATH
        """
        import platform

        if env is None:
            env = os.environ.copy()

        system = platform.system()

        if system == "Darwin":
            if "DYLD_LIBRARY_PATH" not in env:
                ld_path = env.get("LD_LIBRARY_PATH")
                if ld_path:
                    env["DYLD_LIBRARY_PATH"] = ld_path
                    self.logger.debug(
                        "Copied LD_LIBRARY_PATH to DYLD_LIBRARY_PATH for macOS compatibility"
                    )
        elif system == "Linux":
            if "LD_LIBRARY_PATH" not in env:
                dyld_path = env.get("DYLD_LIBRARY_PATH")
                if dyld_path:
                    env["LD_LIBRARY_PATH"] = dyld_path
                    self.logger.debug(
                        "Copied DYLD_LIBRARY_PATH to LD_LIBRARY_PATH for Linux compatibility"
                    )

        return env

    def run(
        self,
        cmd: List[str],
        stdout_callback: Optional[Callable[[str], None]] = None,
        stderr_callback: Optional[Callable[[str], None]] = None,
        expected_returncodes: Optional[set] = None,
        check: bool = False,
    ) -> ProcessResult:
        """
        运行外部命令

        Args:
            cmd: 命令及参数列表
            stdout_callback: 标准输出回调函数
            stderr_callback: 标准错误回调函数
            expected_returncodes: 视为成功的返回码集合，默认仅 {0}。
                某些外部程序（如 PATHd8）会用 returncode=1 表示警告但仍生成输出，
                调用方可传入 {0, 1} 避免日志中出现“command failed”字样。
            check: 返回码不在 ``expected_returncodes`` 中时是否抛出
                ``ExecutionError``（语义同 ``subprocess.run(check=True)``）。

                **默认为 False（只 warning、照常返回 ProcessResult），这是有意
                为之的向后兼容选择**（审阅项 C-35）。理由：现存调用点里有若
                干处按设计需要拿到"非零码 + 那份输出"，直接改默认值会让它们
                从"降级/回报失败对象"变成"抛异常"：

                - ``treepl_method._run_prime``：returncode≠0 时返回
                  ``PrimeResult(success=False)``，由上层继续走默认 opt/optad；
                - ``treepl_method._auto_adjust_opt_params``：CV 探测失败时
                  ``break`` 退出调参循环而非终止分析；
                - ``mcmctree_method`` 速率预估：baseml/codeml 失败时回退到默认
                  ``rgene_gamma`` 先验；
                - ``pathd8_method``：上游成功也可能退出 1，已显式传
                  ``expected_returncodes={0, 1}``；``r8s_pyr8s_method`` 同理
                  存在奇异退出码。

                新写的调用点应当传 ``check=True``（把"必须记得检查 returncode"
                变成"默认安全"）；已在 ``expected_returncodes`` 中列出的码不会
                触发异常，因此容忍特殊码与强制检查可以共存。

        Note:
            子进程的 ``stdin`` 恒接到 ``subprocess.DEVNULL``。r8s 一类的交互式
            shell 在没有 ``-b`` 批处理标志时按手册所述会 "sit there patiently
            waiting for keyboard input"（审阅项 A-6），父进程若不接管 stdin，
            它就一直读到超时为止。给出 EOF 是与命令行标志无关的一层兜底防护：
            即便某个调用点漏了 ``-b``，也只会立刻失败而不是一起挂死。

        Returns:
            ProcessResult: 执行结果

        Raises:
            ExecutionError: 执行失败；或在 ``check=True`` 时返回码不在
                ``expected_returncodes`` 中
            ExecutionTimeoutError: 执行超时
        """
        if expected_returncodes is None:
            expected_returncodes = {0}

        start_time = time.time()
        self._peak_memory = 0.0
        self._stop_monitoring = threading.Event()

        self.logger.info(f"Executing: {' '.join(cmd)}")
        if self.cwd:
            self.logger.info(f"Working directory: {self.cwd}")

        try:
            # text 模式下 bufsize 表示行数，使用默认值 -1（系统默认缓冲）
            self._process = subprocess.Popen(
                cmd,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=-1,
                cwd=self.cwd,
                env=self.env,
            )
            # 嵌套函数 read_stream 里拿不到 ``self._process`` 的收窄（它可被其他
            # 方法置 None），所以先取一个本地引用：下面的终止/取流/等待都用它。
            process = self._process

            if self.monitor_memory:
                self._start_memory_monitor(self._process.pid)

            stdout_lines: List[str] = []
            stderr_lines: List[str] = []
            # C-34：计数与超限判定在锁内合成一步，避免两个读线程互相覆盖增量
            output_budget = _OutputBudget(self.max_output_size)

            def read_stream(
                stream: IO[str],
                callback: Optional[Callable[[str], None]],
                storage: List[str],
                is_stdout: bool = True,
            ) -> None:
                for line in iter(stream.readline, ""):
                    if not line:
                        break
                    line = line.rstrip("\n")

                    if output_budget.record(len(line) + 1):
                        self.logger.error(
                            f"Output size exceeded maximum limit "
                            f"({self.max_output_size} bytes)"
                        )
                        process.terminate()
                        return

                    storage.append(line)
                    if callback:
                        callback(line)
                    else:
                        if is_stdout:
                            self.logger.info(line)
                        else:
                            self.logger.error(line)
                stream.close()

            stdout_thread = threading.Thread(
                target=read_stream,
                args=(process.stdout, stdout_callback, stdout_lines, True),
            )
            stderr_thread = threading.Thread(
                target=read_stream,
                args=(process.stderr, stderr_callback, stderr_lines, False),
            )

            stdout_thread.start()
            stderr_thread.start()

            self._register_signal_handlers()

            try:
                returncode = process.wait(timeout=self.timeout)
            except subprocess.TimeoutExpired:
                self._log_intermediate_results()
                self._cleanup()
                stdout_thread.join(timeout=5)
                stderr_thread.join(timeout=5)
                raise ExecutionTimeoutError(
                    f"Command timed out after {self.timeout}s: {' '.join(cmd)}. "
                    f"Intermediate results (if any) have been logged. "
                    f"Possible causes: parameter error causing infinite loop, or insufficient burn-in."
                )
            finally:
                self._unregister_signal_handlers()

            stdout_thread.join()
            stderr_thread.join()

            # 检查 reader 线程是否因输出超限而终止
            if output_budget.exceeded:
                self._cleanup()
                raise output_budget.error or ExecutionError(
                    f"Output size exceeded maximum limit ({self.max_output_size} bytes)"
                )

            self._cleanup()

            execution_time = time.time() - start_time

            stdout = "\n".join(stdout_lines)
            stderr = "\n".join(stderr_lines)

            self.logger.info(
                f"Process completed in {execution_time:.2f}s "
                f"(returncode={returncode})"
            )

            if returncode not in expected_returncodes:
                from ..core.exceptions import classify_execution_error

                category, suggestion = classify_execution_error(stderr, returncode)
                self.logger.warning(
                    f"External command failed with returncode {returncode}. "
                    f"Category: {category.value}. Suggestion: {suggestion}"
                )
                if check:
                    # C-35：默认安全的失败传播。异常消息保留分类结论与 stderr 末尾，
                    # 使调用方即使不读日志也能拿到诊断依据。
                    stderr_tail = "\n".join(stderr.splitlines()[-10:])
                    message = (
                        f"Command exited with unexpected returncode {returncode} "
                        f"(expected one of {sorted(expected_returncodes)}): "
                        f"{' '.join(cmd)}\n"
                        f"Category: {category.value}. Suggestion: {suggestion}"
                    )
                    if stderr_tail:
                        message += f"\n--- last stderr lines ---\n{stderr_tail}"
                    raise ExecutionError(message)
            elif returncode != 0:
                self.logger.info(
                    f"External command returned {returncode} (within expected range); "
                    f"continuing with output parsing."
                )

            if self._peak_memory > 0:
                self.logger.info(f"Peak memory: {self._peak_memory:.1f} MB")

            return ProcessResult(
                returncode=returncode,
                stdout=stdout,
                stderr=stderr,
                execution_time=execution_time,
                peak_memory_mb=self._peak_memory if self._peak_memory > 0 else None,
            )

        except ExecutionTimeoutError:
            raise
        except ExecutionError:
            raise
        except MemoryError:
            self._cleanup()
            raise ExecutionError(
                "System memory exhausted. The analysis requires more memory than available.\n"
                "Suggestions:\n"
                "  1. Reduce the number of parallel threads (--threads)\n"
                "  2. Use a machine with more RAM\n"
                "  3. Reduce dataset size or complexity\n"
                "  4. For MCMC methods, reduce nsample or burnin values"
            )
        except Exception as e:
            self._cleanup()
            error_str = str(e).lower()
            if any(
                keyword in error_str
                for keyword in ["memory", "oom", "killed", "signal 9"]
            ):
                raise ExecutionError(
                    f"Process terminated, likely due to memory exhaustion: {e}\n"
                    "Suggestions:\n"
                    "  1. Reduce the number of parallel threads (--threads)\n"
                    "  2. Use a machine with more RAM\n"
                    "  3. For large datasets, consider using methods with lower "
                    "memory requirements (e.g., PATHd8 instead of MCMCTree)"
                ) from e
            raise ExecutionError(f"Failed to execute command: {e}") from e

    def _register_signal_handlers(self) -> None:
        """注册信号处理器。

        注意：signal.signal() 只能在主线程调用。当 ProcessRunner 运行在
        线程池（futures.thread）中时，跳过信号注册，避免
        `signal only works in main thread of the main interpreter`。
        """
        if threading.current_thread() is not threading.main_thread():
            self.logger.debug(
                "Skipping signal-handler registration: not running in main thread "
                f"({threading.current_thread().name})"
            )
            self._original_sigterm_handler = None
            self._original_sigint_handler = None
            return

        def signal_handler(signum: int, frame: Optional[FrameType]) -> None:
            self.logger.warning(f"Received signal {signum}, terminating process")
            self._cleanup()
            if signum == signal.SIGINT:
                # 用户中断需要以 KeyboardInterrupt 传播，以便上层返回 130
                raise KeyboardInterrupt("Process interrupted by user")
            raise ExecutionError(f"Process interrupted by signal {signum}")

        self._original_sigterm_handler = signal.signal(signal.SIGTERM, signal_handler)
        self._original_sigint_handler = signal.signal(signal.SIGINT, signal_handler)

    def _unregister_signal_handlers(self) -> None:
        """恢复原始信号处理器"""
        if self._original_sigterm_handler is not None:
            signal.signal(signal.SIGTERM, self._original_sigterm_handler)
            self._original_sigterm_handler = None
        if self._original_sigint_handler is not None:
            signal.signal(signal.SIGINT, self._original_sigint_handler)
            self._original_sigint_handler = None

    def _cleanup(self) -> None:
        """清理资源（可重入安全）"""
        if self._stop_monitoring:
            self._stop_monitoring.set()

        process = self._process
        self._process = None
        if process:
            try:
                process.terminate()
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                try:
                    process.kill()
                except Exception:
                    pass
            except Exception:
                pass

        if self._memory_monitor_thread and self._memory_monitor_thread.is_alive():
            self._memory_monitor_thread.join(timeout=3)
            if self._memory_monitor_thread.is_alive():
                self.logger.warning("Memory monitor thread did not terminate in time")

    def _log_intermediate_results(self) -> None:
        """
        在超时发生时记录中间结果

        输出当前已完成的检查点信息，帮助诊断超时原因
        """
        if self.cwd and self.cwd.exists():
            intermediate_files = list(self.cwd.glob("*"))
            self.logger.warning(
                f"Timeout occurred. Current working directory: {self.cwd}. "
                f"Intermediate files found: {len(intermediate_files)}"
            )
            for f in intermediate_files[:20]:
                size = f.stat().st_size if f.is_file() else 0
                self.logger.warning(f"  - {f.name} ({size} bytes)")
            if len(intermediate_files) > 20:
                self.logger.warning(
                    f"  ... and {len(intermediate_files) - 20} more files"
                )

    def _start_memory_monitor(self, pid: int) -> None:
        """启动内存监控线程"""
        try:
            import psutil
        except ImportError:
            self.logger.warning("psutil not available, memory monitoring disabled")
            return

        # 停止事件由 run() 在调本方法前建好。没有它就没法通知监控线程退出，
        # 必须当场报错，而不是悄悄起一个停不下来的线程。
        stop_event = self._stop_monitoring
        if stop_event is None:
            raise RuntimeError(
                "ProcessRunner._start_memory_monitor() called before run() created "
                "the stop event; the monitor thread would be unstoppable"
            )
        stop_event.clear()

        def monitor() -> None:
            try:
                proc = psutil.Process(pid)
                while not stop_event.is_set():
                    try:
                        mem_mb = proc.memory_info().rss / 1024 / 1024
                        self._peak_memory = max(self._peak_memory, mem_mb)

                        total_mem = psutil.virtual_memory().total / 1024 / 1024
                        available_mem = psutil.virtual_memory().available / 1024 / 1024

                        if mem_mb > total_mem * 0.95:
                            self.logger.error(
                                f"CRITICAL: Memory usage ({mem_mb:.1f}MB) exceeds 95% of system memory! "
                                f"Available: {available_mem:.1f}MB. "
                                "Process may be terminated by system OOM killer. "
                                "Immediately reducing workload is recommended."
                            )
                        elif mem_mb > total_mem * 0.8:
                            self.logger.warning(
                                f"Memory usage ({mem_mb:.1f}MB) exceeds 80% of system memory. "
                                f"Available: {available_mem:.1f}MB. "
                                "Consider reducing --threads or using a machine with more RAM."
                            )

                        time.sleep(0.5)
                    except psutil.NoSuchProcess:
                        break
                    except Exception:
                        break
            except Exception as e:
                self.logger.debug(f"Memory monitor error: {e}")

        self._memory_monitor_thread = threading.Thread(target=monitor, daemon=True)
        self._memory_monitor_thread.start()

    def check_executable(self, executable: str) -> bool:
        """检查可执行文件是否可用（跨平台）"""
        return shutil.which(executable) is not None

    def get_version(
        self, executable: str, version_flag: str = "--version"
    ) -> Optional[str]:
        """获取软件版本信息

        优先在 returncode==0 时返回 stdout/stderr 内容；若 returncode 非 0 但
        输出看起来像版本号(例如 mcmctree 在打印版本 banner 后尝试把 ``--version``
        当作输入文件打开而退出非零,但 stdout 含 "paml version 4.10.8"),
        也回退返回最佳可用文本,避免"unknown version"。
        """
        try:
            result = subprocess.run(
                [executable, version_flag],
                capture_output=True,
                text=True,
                timeout=10,
                # 同 run()：不给交互式程序一个"等待键盘输入"的机会（A-6）
                stdin=subprocess.DEVNULL,
            )
            text = (result.stdout or "").strip() or (result.stderr or "").strip()
            if not text:
                return None
            # 取第一个非空行,捕获 "paml version X" / "v1.2" 等
            first_line = next(
                (ln.strip() for ln in text.splitlines() if ln.strip()), text
            )
            if result.returncode == 0:
                return first_line
            # returncode 非 0:仅当首行看起来像版本号时才回退(避开 usage/错误信息)
            import re as _re

            if _re.search(r"version|v?\d+\.\d+", first_line, _re.IGNORECASE):
                return first_line
            return None
        except Exception:
            return None
