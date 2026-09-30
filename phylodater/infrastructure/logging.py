"""
PhyloDaterLogger - 实时流式日志系统

特性：
- 实时输出：使用 sys.stdout 实时刷新，不缓冲
- ISO 8601 时间戳：2025-03-21T10:15:30.123
- 日志格式：[TIMESTAMP] [LEVEL] [MODULE/FUNCTION] message
- 多级别支持：DEBUG、INFO、WARNING、ERROR、CRITICAL
- 禁止彩色输出，纯文本格式
- 日志文件持久化：--log-file 参数，UTF-8 编码，实时 flush
- 路径脱敏：自动将日志中的绝对路径转换为相对路径或掩码
"""

import inspect
import re
import sys
import threading
import time
from datetime import datetime
from enum import Enum
from pathlib import Path
from types import TracebackType
from typing import IO, Any, Dict, List, Literal, Optional, Union

# 调用栈跳过的帧数：_get_caller_module -> _log -> debug/info/... -> 实际调用者。
# 提取为模块级常量，便于按需调整而非硬编码在循环里。
CALLER_FRAME_SKIP = 4


class LogLevel(Enum):
    """日志级别"""

    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class StepTimer:
    """步骤计时器"""

    def __init__(self) -> None:
        self._start_time = time.time()
        self._step_start = time.time()
        self._steps: List[Dict[str, Any]] = []

    def start_step(self, step_name: str) -> None:
        """开始一个新步骤"""
        self._step_start = time.time()

    def end_step(self, step_name: str) -> float:
        """结束当前步骤，返回耗时"""
        elapsed = time.time() - self._step_start
        self._steps.append(
            {
                "name": step_name,
                "duration": elapsed,
                "timestamp": datetime.now().isoformat(timespec="milliseconds"),
            }
        )
        return elapsed

    def get_total_elapsed(self) -> float:
        """获取总耗时"""
        return time.time() - self._start_time

    def get_step_summary(self) -> List[Dict[str, Any]]:
        """获取步骤摘要"""
        return self._steps.copy()


class PhyloDaterLogger:
    """
    实时流式日志管理器

    特性：
    - 实时输出：使用 sys.stdout 实时刷新，不缓冲
    - ISO 8601 时间戳格式
    - 模块/函数名追踪
    - 多级别支持：DEBUG、INFO、WARNING、ERROR、CRITICAL
    - 禁止彩色输出
    """

    _instance = None
    _lock = threading.Lock()
    # 只标注、不给类级默认值：两个字段仍然由 __new__ / __init__ 在实例上写入，
    # 语义零变化。不标注时 mypy 推不出类型（has-type），而这两个属性
    # 是单例防重入与步骤计时器的核心状态，必须真的被检查一下。
    _initialized: bool
    _step_start_time: Optional[float]

    def __new__(cls, *args: Any, **kwargs: Any) -> "PhyloDaterLogger":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._initialized = False
        return cls._instance

    def __init__(
        self, log_file: Optional[Path] = None, level: Optional[LogLevel] = None
    ) -> None:
        if self._initialized:
            # 单例已存在时，允许更新日志文件路径和级别（例如 CLI 解析出 --log-file 后）
            if log_file is not None:
                self._update_log_file(log_file)
            if level is not None:
                self.level = level
            return

        self.log_file = Path(log_file) if log_file is not None else None
        self.level = level if level is not None else LogLevel.INFO
        self._base_dir = Path.cwd()
        self._timer = StepTimer()
        self._current_step = 0
        self._step_start_time = None

        # 敏感路径模式
        self._sensitive_patterns = [
            re.compile(r"/home/[^/]+/"),
            re.compile(r"/Users/[^/]+/"),
            re.compile(r"/tmp/[^/]+/"),
            re.compile(r"/var/[^/]+/"),
        ]

        # 日志文件句柄（未启用文件日志时为 None）
        self._log_file_handle: Optional[IO[str]] = None

        self._initialized = True

    def _update_log_file(self, log_file: Union[str, Path]) -> None:
        """在单例已初始化后更新日志文件路径（用于 CLI 晚于模块导入指定 --log-file）。"""
        new_log_file = Path(log_file)
        if (
            self.log_file is not None
            and self.log_file.resolve() == new_log_file.resolve()
        ):
            return
        self.log_file = new_log_file
        # 如果已经 start() 过，需要关闭旧句柄并打开新文件
        if self._log_file_handle:
            try:
                self._log_file_handle.flush()
                self._log_file_handle.close()
            except Exception:
                pass
            self._log_file_handle = None
            self.start()

    def start(self) -> None:
        """启动日志系统"""
        if self._log_file_handle:
            return

        # 设置 stdout 为无缓冲模式
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(write_through=True)

        # 创建日志文件
        if self.log_file:
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            self._log_file_handle = open(
                self.log_file, "a", encoding="utf-8", newline="\n"
            )

        self._log(LogLevel.INFO, "日志系统已启动", module="logging")

    def stop(self) -> None:
        """停止日志系统"""
        if self._log_file_handle:
            self._log_file_handle.flush()
            self._log_file_handle.close()
            self._log_file_handle = None

    def set_base_dir(self, base_dir: Path) -> None:
        """设置基准目录用于相对路径计算"""
        self._base_dir = base_dir

    def _sanitize_path(self, path: str) -> str:
        """脱敏单个路径（任何失败都原样返回，日志层绝不因脱敏而抛异常）。"""
        try:
            path_obj = Path(path)

            try:
                rel_path = path_obj.resolve().relative_to(self._base_dir.resolve())
                # 当路径就是基准目录本身时，保留原始路径而不是显示 "."
                if rel_path == Path("."):
                    return path
                return str(rel_path)
            except ValueError:
                pass

            for pattern in self._sensitive_patterns:
                if pattern.search(path):
                    return pattern.sub("***/", path)

            return path
        except Exception:
            return path

    @staticmethod
    def _as_text(message: Any) -> str:
        """把任意日志载荷转成字符串。

        审阅项 B-27 [已证]：``_sanitize_message`` 对消息做 ``pattern.sub(repl, msg)``，
        ``msg`` 必须是字符串。适配器里 ``logger.warning(some_exception_object)``
        （"报告降级的代码"）会直接抛 ``TypeError: expected string or bytes-like
        object`` ——记录问题的通道自己制造了新问题，而且错误文案与真实成因相距极远。
        这里在**框架层**兜底：任何非字符串（Exception、Warning、Path、dict…）先
        ``str()`` 再进正则。
        """
        if isinstance(message, str):
            return message
        if message is None:
            return "None"
        try:
            return str(message)
        except Exception:  # 自定义 __str__ 本身崩溃时的最后退路（日志层绝不抛）
            return f"<unprintable {type(message).__name__} object>"

    def _sanitize_message(self, msg: Any) -> str:
        """脱敏消息中的路径（非字符串消息先转成字符串，见 ``_as_text``）。"""
        msg = self._as_text(msg)
        path_pattern = re.compile(r"(?:[A-Za-z]:\\|/)(?:[^/\s]+[/\\])+(?:[^/\s]+)")

        def replace_path(match: "re.Match[str]") -> str:
            return self._sanitize_path(match.group(0))

        return path_pattern.sub(replace_path, msg)

    def _format_timestamp(self) -> str:
        """
        格式化时间戳（ISO 8601）

        格式: 2025-03-21T10:15:30.123
        """
        return datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]

    def _format_duration(self, seconds: float) -> str:
        """格式化持续时间"""
        if seconds < 0.001:
            return "< 0.001s"
        elif seconds < 1:
            return f"{seconds*1000:.1f}ms"
        elif seconds < 60:
            return f"{seconds:.2f}s"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            return f"{minutes}m {secs:.1f}s"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            secs = seconds % 60
            return f"{hours}h {minutes}m {secs:.1f}s"

    def _get_caller_module(self) -> str:
        """获取调用者的模块/函数名"""
        try:
            # 向上查找调用栈，跳过日志系统自身的帧
            frame = inspect.currentframe()
            # 跳过: _get_caller_module -> _log -> debug/info/... -> 实际调用者
            for _ in range(CALLER_FRAME_SKIP):
                if frame and frame.f_back:
                    frame = frame.f_back
                else:
                    break

            if frame:
                module = frame.f_globals.get("__name__", "?")
                func = frame.f_code.co_name
                # 简化模块名：只取最后两级
                parts = module.split(".")
                if len(parts) > 2:
                    module = ".".join(parts[-2:])
                return f"{module}/{func}"
        except Exception:
            pass
        return "?"

    def _output(self, message: str, is_file_only: bool = False) -> None:
        """输出消息到控制台和文件"""
        # 输出到控制台
        if not is_file_only:
            print(message, flush=True)

        # 输出到文件（同样格式，实时flush）
        if self._log_file_handle:
            self._log_file_handle.write(message + "\n")
            self._log_file_handle.flush()

    def _log(
        self,
        level: LogLevel,
        message: Any,
        is_header: bool = False,
        indent: int = 0,
        is_step: bool = False,
        is_kv: bool = False,
        module: Optional[str] = None,
        is_file_only: bool = False,
    ) -> None:
        """内部日志方法。

        ``message`` 声明为 ``Any``：所有公共入口（debug/info/warning/...）都接受
        任意对象并在此统一转字符串，日志层本身绝不因载荷类型而抛异常（B-27）。
        """
        message = self._as_text(message)
        # 检查日志级别
        level_order = [
            LogLevel.DEBUG,
            LogLevel.INFO,
            LogLevel.WARNING,
            LogLevel.ERROR,
            LogLevel.CRITICAL,
        ]
        if level_order.index(level) < level_order.index(self.level):
            return

        # 脱敏消息（非字符串载荷在这里被安全地转成字符串）
        sanitized_msg = self._sanitize_message(message)

        # 格式化时间戳 (ISO 8601)
        timestamp = self._format_timestamp()

        # 获取级别名称（右对齐9字符）
        level_str = level.value.ljust(9)

        # 获取模块/函数名
        if module is None:
            module = self._get_caller_module()

        # 构建消息
        if is_header:
            # 章节标题
            line = f"\n{'='*70}"
            self._output(line)
            line = f"  {sanitized_msg}"
            self._output(line)
            line = f"{'='*70}"
            self._output(line)
        elif is_step:
            # 步骤标题
            self._current_step += 1
            self._step_start_time = time.time()
            step_num = self._current_step
            line = f"\n[STEP {step_num}] {sanitized_msg}"
            self._output(line)
            line = f"{'─'*70}"
            self._output(line)
        elif is_kv:
            # 键值对
            line = f"{timestamp} | {level_str} |     {sanitized_msg}"
            self._output(line)
        else:
            # 普通消息
            # 格式: 2025-03-21T10:15:30.123 | INFO     | module/func message
            indent_str = "  " * indent
            line = f"{timestamp} | {level_str} | {module} {sanitized_msg}"
            self._output(f"{indent_str}{line}", is_file_only=is_file_only)

    # ==================== 公共接口 ====================
    #
    # 以下所有 ``message`` 都接受任意对象（Exception / Warning / Path / dict…）：
    # 载荷在进入正则脱敏前统一 ``str()``，日志层绝不因消息类型抛异常（B-27）。

    def debug(self, message: Any, **kwargs: Any) -> None:
        """记录调试信息"""
        module = kwargs.pop("module", None)
        self._log(LogLevel.DEBUG, message, module=module)

    def debug_file_only(self, message: Any) -> None:
        """仅写入日志文件的调试信息（不输出到控制台）"""
        self._log(LogLevel.DEBUG, message, is_file_only=True)

    def info(self, message: Any, **kwargs: Any) -> None:
        """记录常规信息"""
        module = kwargs.pop("module", None)
        self._log(LogLevel.INFO, message, module=module)

    def warning(self, message: Any, **kwargs: Any) -> None:
        """记录警告信息（可选地同时通过 warnings 发出一个类别警告）"""
        import warnings as _warnings

        module = kwargs.pop("module", None)
        category = kwargs.pop("category", None)
        text = self._as_text(message)
        if category is not None:
            try:
                _warnings.warn(text, category, stacklevel=2)
            except Exception as exc:  # 发警告的通道也不能自己炸掉
                self._log(
                    LogLevel.WARNING,
                    f"{text} (warnings.warn({category!r}) failed: {exc})",
                    module=module,
                )
                return
        self._log(LogLevel.WARNING, text, module=module)

    def error(self, message: Any, **kwargs: Any) -> None:
        """记录错误信息"""
        module = kwargs.pop("module", None)
        self._log(LogLevel.ERROR, message, module=module)

    def critical(self, message: Any, **kwargs: Any) -> None:
        """记录严重错误"""
        module = kwargs.pop("module", None)
        self._log(LogLevel.CRITICAL, message, module=module)

    def section(self, title: Any) -> None:
        """打印章节标题"""
        self._log(LogLevel.INFO, title, is_header=True, module="main")

    def step(self, step_num: int, title: Any) -> None:
        """打印步骤标题"""
        self._current_step = step_num - 1  # 会被 _log 中的 is_step 递增
        self._log(LogLevel.INFO, title, is_step=True, module="main")

    def kv(self, key: Any, value: Any, indent: int = 4) -> None:
        """打印键值对"""
        self._log(
            LogLevel.INFO,
            f"{key!s:<25}: {value}",
            is_kv=True,
            indent=1,
            module="config",
        )

    def success(self, message: Any) -> None:
        """打印成功消息"""
        timestamp = self._format_timestamp()
        line = f"{timestamp} | {'PASS'.ljust(9)} | {self._as_text(message)}"
        self._output(line)

    def fail(self, message: Any) -> None:
        """打印失败消息"""
        timestamp = self._format_timestamp()
        line = f"{timestamp} | {'FAIL'.ljust(9)} | {self._as_text(message)}"
        self._output(line)

    def progress(self, current: int, total: int, message: str = "") -> None:
        """打印进度"""
        if total > 0:
            percent = (current / total) * 100
            bar_length = 30
            filled = int(bar_length * current / total)
            bar = "=" * filled + "-" * (bar_length - filled)

            line = f"\r  [{bar}] {percent:5.1f}% ({current}/{total})"

            if message:
                line += f" {message}"

            # 不换行输出
            print(line, end="", flush=True)

            if current >= total:
                print()  # 完成时换行

    def timing(self, operation: str, duration: float) -> None:
        """打印计时信息"""
        duration_str = self._format_duration(duration)
        timestamp = self._format_timestamp()
        line = f"{timestamp} | {'TIMER'.ljust(9)} | {operation}: {duration_str}"
        self._output(line)

    def step_complete(self, message: str = "") -> None:
        """标记当前步骤完成"""
        if self._step_start_time:
            elapsed = time.time() - self._step_start_time
            duration_str = self._format_duration(elapsed)
            timestamp = self._format_timestamp()
            line = f"{timestamp} | {'PASS'.ljust(9)} | Step complete ({duration_str})"
            if message:
                line += f" - {message}"
            self._output(line)
            self._step_start_time = None

    def execution_provenance(
        self,
        tool_name: str,
        command: str,
        version: str,
        environment: Dict[str, Any],
        artifacts: Dict[str, Any],
    ) -> None:
        """记录工具执行的血统信息"""
        self.info(f"执行 {tool_name} v{version}", module="provenance")
        self.debug(f"命令: {command}", module="provenance")
        if artifacts:
            for key, value in artifacts.items():
                self.debug(f"  {key}: {value}", module="provenance")

    def __enter__(self) -> "PhyloDaterLogger":
        self.start()
        return self

    def __exit__(
        self,
        exc_type: Optional[type],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> Literal[False]:
        self.stop()
        return False


def get_logger(
    log_file: Optional[Path] = None, level: Optional[LogLevel] = None
) -> PhyloDaterLogger:
    """获取日志记录器实例。

    注意：单例的日志文件路径和级别可以通过本函数更新；不传参数时保持现有
    配置不变，避免覆盖 CLI 已设置的 --log-file/--log-level。
    """
    return PhyloDaterLogger(log_file, level)


def setup_logging(
    level: LogLevel = LogLevel.INFO, log_file: Optional[Path] = None
) -> PhyloDaterLogger:
    """设置日志系统"""
    logger = PhyloDaterLogger(log_file, level)
    logger.start()
    return logger
