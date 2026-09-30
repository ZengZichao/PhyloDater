"""
SignalHandler - 信号处理与优雅终止

注册 SIGINT (Ctrl+C) 和 SIGTERM 信号处理器，实现：
- 优雅终止：停止当前操作，输出进度，清理资源
- 临时文件清理：删除未完成的临时输出
- 资源释放：确保日志文件句柄正确关闭
"""

import atexit
import os
import shutil
import signal
import sys
import threading
from dataclasses import dataclass, field
from pathlib import Path
from types import FrameType, TracebackType
from typing import Any, Callable, Literal, Optional, Set

from .logging import get_logger


@dataclass
class ShutdownState:
    """关闭状态跟踪"""

    is_shutting_down: bool = False
    shutdown_reason: Optional[str] = None
    progress_info: dict = field(default_factory=dict)
    temp_dirs: Set[Path] = field(default_factory=set)
    temp_files: Set[Path] = field(default_factory=set)
    cleanup_callbacks: list = field(default_factory=list)
    original_sigint: Any = None
    original_sigterm: Any = None
    #: 信号是否已注册。旧代码只在 ``register()`` 里随手 ``self._state.registered
    #: = True``，字段根本没声明，读侧只能靠 ``getattr(..., False)`` 兑；
    #: 把它补成字段才是这个状态机的真实形状。
    registered: bool = False
    checkpoint_manager: Optional[Any] = None


class GracefulShutdownHandler:
    """
    优雅关闭处理器

    注册信号处理器，实现优雅终止和资源清理。
    """

    # 单例模式
    _instance: Optional["GracefulShutdownHandler"] = None
    _lock = threading.Lock()
    # 只标注不给类级默认值：仍由 __new__ / __init__ 在实例上写。
    _initialized: bool

    def __new__(cls) -> "GracefulShutdownHandler":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._initialized = False
        return cls._instance

    def __init__(self) -> None:
        if self._initialized:
            return

        self._state = ShutdownState()
        self._logger = get_logger()
        self._initialized = True

    @property
    def is_shutting_down(self) -> bool:
        """是否正在关闭"""
        return self._state.is_shutting_down

    @property
    def progress_info(self) -> dict:
        """获取进度信息"""
        return self._state.progress_info.copy()

    def register(self) -> None:
        """注册信号处理器"""
        if getattr(self._state, "registered", False):
            self._logger.debug("Signal handlers already registered; skipping")
            return

        # 保存原始信号处理器
        self._state.original_sigint = signal.getsignal(signal.SIGINT)
        self._state.original_sigterm = signal.getsignal(signal.SIGTERM)

        # 注册新的信号处理器
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

        # 注册退出清理
        atexit.register(self._cleanup_on_exit)

        self._state.registered = True
        self._logger.debug("Signal handlers registered (SIGINT, SIGTERM)")

    def unregister(self) -> None:
        """恢复原始信号处理器"""
        if not self._state.registered:
            return

        if self._state.original_sigint is not None:
            signal.signal(signal.SIGINT, self._state.original_sigint)
        if self._state.original_sigterm is not None:
            signal.signal(signal.SIGTERM, self._state.original_sigterm)

        # 注销退出清理
        try:
            atexit.unregister(self._cleanup_on_exit)
        except AttributeError:
            pass

        self._state.registered = False
        self._state.original_sigint = None
        self._state.original_sigterm = None
        self._logger.debug("Signal handlers unregistered")

    def _signal_handler(self, signum: int, frame: Optional[FrameType]) -> None:
        """信号处理函数"""
        if self._state.is_shutting_down:
            # 如果已经在关闭过程中，强制退出
            self._logger.warning("Forced exit requested during shutdown")
            sys.exit(130)

        self._state.is_shutting_down = True

        # 确定信号类型
        if signum == signal.SIGINT:
            self._state.shutdown_reason = "User interrupt (Ctrl+C)"
            sig_name = "SIGINT"
        elif signum == signal.SIGTERM:
            self._state.shutdown_reason = "Termination signal"
            sig_name = "SIGTERM"
        else:
            self._state.shutdown_reason = f"Signal {signum}"
            sig_name = f"Signal {signum}"

        self._logger.warning(f"\n{'='*60}")
        self._logger.warning(f"Received {sig_name} - Initiating graceful shutdown...")
        self._logger.warning(f"{'='*60}")

        # 输出当前进度
        self._print_progress()

        # 保存检查点进度
        checkpoint_mgr = getattr(self._state, "checkpoint_manager", None)
        if checkpoint_mgr is not None:
            try:
                checkpoint_mgr.save()
                self._logger.info("Checkpoint progress saved.")
            except Exception as e:
                self._logger.warning(f"Failed to save checkpoint: {e}")

        # 执行清理
        self._cleanup()

        # 退出
        self._logger.warning("Shutdown complete.")
        sys.exit(130)

    def _print_progress(self) -> None:
        """输出当前进度"""
        if not self._state.progress_info:
            self._logger.info("No progress information available.")
            return

        self._logger.info("Current progress:")

        info = self._state.progress_info

        if "current_step" in info:
            self._logger.info(f"  Current step: {info['current_step']}")

        if "trees_processed" in info:
            self._logger.info(f"  Trees processed: {info['trees_processed']}")

        if "trees_total" in info:
            self._logger.info(f"  Trees total: {info['trees_total']}")

        if "methods_completed" in info:
            self._logger.info(f"  Methods completed: {info['methods_completed']}")

        if "methods_total" in info:
            self._logger.info(f"  Methods total: {info['methods_total']}")

        if "elapsed_time" in info:
            elapsed = info["elapsed_time"]
            if elapsed < 60:
                self._logger.info(f"  Elapsed time: {elapsed:.1f}s")
            elif elapsed < 3600:
                self._logger.info(f"  Elapsed time: {elapsed/60:.1f}m")
            else:
                self._logger.info(f"  Elapsed time: {elapsed/3600:.1f}h")

    def _cleanup(self) -> None:
        """执行清理操作"""
        self._logger.info("Cleaning up resources...")

        # 删除临时文件
        for temp_file in self._state.temp_files:
            try:
                if temp_file.exists():
                    temp_file.unlink()
                    self._logger.debug(f"Removed temporary file: {temp_file}")
            except Exception as e:
                self._logger.warning(
                    f"Failed to remove temporary file {temp_file}: {e}"
                )

        # 删除临时目录
        for temp_dir in self._state.temp_dirs:
            try:
                if temp_dir.exists():
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    self._logger.debug(f"Removed temporary directory: {temp_dir}")
            except Exception as e:
                self._logger.warning(
                    f"Failed to remove temporary directory {temp_dir}: {e}"
                )

        # 执行注册的清理回调
        for callback in self._state.cleanup_callbacks:
            try:
                callback()
            except Exception as e:
                self._logger.warning(f"Cleanup callback failed: {e}")

    def _cleanup_on_exit(self) -> None:
        """退出时的清理（atexit回调）"""
        if not self._state.is_shutting_down:
            # 正常退出时也执行清理
            self._cleanup()

    def update_progress(self, **kwargs: Any) -> None:
        """更新进度信息"""
        self._state.progress_info.update(kwargs)

    def register_temp_dir(self, path: Path) -> None:
        """注册临时目录用于清理"""
        self._state.temp_dirs.add(path)
        # 显式设置临时目录权限为 0o700，确保仅当前用户可访问
        try:
            os.chmod(str(path), 0o700)
        except Exception as e:
            self._logger.warning(f"Failed to chmod temp directory {path}: {e}")

    def unregister_temp_dir(self, path: Path) -> None:
        """取消注册临时目录"""
        self._state.temp_dirs.discard(path)

    def register_temp_file(self, path: Path) -> None:
        """注册临时文件用于清理"""
        self._state.temp_files.add(path)

    def unregister_temp_file(self, path: Path) -> None:
        """取消注册临时文件"""
        self._state.temp_files.discard(path)

    def register_cleanup_callback(self, callback: Callable) -> None:
        """注册清理回调函数"""
        self._state.cleanup_callbacks.append(callback)

    def unregister_cleanup_callback(self, callback: Callable) -> None:
        """取消注册清理回调函数"""
        try:
            self._state.cleanup_callbacks.remove(callback)
        except ValueError:
            pass


def get_shutdown_handler() -> GracefulShutdownHandler:
    """获取关闭处理器单例"""
    return GracefulShutdownHandler()


class ShutdownContext:
    """
    关闭上下文管理器

    用于在 with 语句中自动注册和注销信号处理器。
    """

    def __init__(self, output_dir: Optional[Path] = None) -> None:
        self._handler = get_shutdown_handler()
        self._output_dir = output_dir
        # 独立的临时目录（由 __enter__ 创建），绝不指向用户的 output_dir，
        # 以避免在清理时递归删除整个结果目录。
        self._temp_dir: Optional[Path] = None

    def __enter__(self) -> "GracefulShutdownHandler":
        self._handler.register()

        # 注册独立的临时目录用于清理。
        # 关键修复：绝不能把用户的 output_dir 注册为待清理目录——否则在 Ctrl+C
        # 或 atexit 时会通过 shutil.rmtree 递归删除整个结果目录，造成数据丢失。
        # 改为在系统临时区创建一个独立的临时目录，仅清理该子目录。
        if self._output_dir:
            import tempfile

            self._temp_dir = Path(tempfile.mkdtemp(prefix="phylodater_"))
            self._handler.register_temp_dir(self._temp_dir)

        return self._handler

    def __exit__(
        self,
        exc_type: Optional[type],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> Literal[False]:
        self._handler.unregister()
        return False
