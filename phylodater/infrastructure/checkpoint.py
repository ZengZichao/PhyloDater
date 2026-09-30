"""
Checkpoint Manager - 检查点/恢复功能

支持长时间运行任务（如 MCMCTree）的断点续传
"""

import json
import os
import threading
from dataclasses import asdict, dataclass, field
from datetime import datetime
from enum import Enum, auto
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Mapping, Optional, Union

if TYPE_CHECKING:
    from ..models import DatingResult

from .logging import get_logger
from .safe_io import safe_writer

#: 运行级"输入指纹"在 ``input_hashes`` 里占用的逻辑键（不是真实路径）。
#: ``core/pipeline.py`` 的 B-23 接线用的就是同名常量；此处独立定义以避免
#: infrastructure → core 的反向依赖，两者必须保持一致。
RUN_FINGERPRINT_KEY = "phylodater:run-fingerprint"

#: 读取"上次运行的输入指纹"时依次尝试的逻辑键（兼容不同的登记口径）。
_FINGERPRINT_KEYS = (
    RUN_FINGERPRINT_KEY,
    "run_fingerprint",
    "__run_fingerprint__",
    "__aggregate__",
)


class CheckpointStatus(Enum):
    """检查点状态"""

    PENDING = auto()  # 等待执行
    RUNNING = auto()  # 正在执行
    COMPLETED = auto()  # 已完成
    FAILED = auto()  # 失败
    SKIPPED = auto()  # 已跳过


@dataclass
class StepCheckpoint:
    """单个步骤的检查点信息"""

    step_name: str
    status: str
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    #: 默认空字典（``default_factory``），而不是把 ``None`` 标成 ``Dict`` 再靠
    #: ``__post_init__`` 补平：后者在静态视角下会让任何给 ``data`` 的赋值成类型错。
    data: Dict[str, Any] = field(default_factory=dict)
    error_message: Optional[str] = None

    def __post_init__(self) -> None:
        if self.data is None:
            self.data = {}


class CheckpointManager:
    """
    检查点管理器

    功能：
    - 保存和加载检查点
    - 支持步骤级别的恢复
    - 自动清理过期检查点
    - 线程安全
    """

    def __init__(
        self, base_dir: Union[str, Path], pipeline_name: str = "pipeline"
    ) -> None:
        """
        初始化检查点管理器

        Args:
            base_dir: 基础目录路径（支持字符串或 Path）
            pipeline_name: 流程名称
        """
        # 支持字符串或 Path 类型
        self.base_dir = Path(base_dir) if isinstance(base_dir, str) else base_dir

        # 路径遍历防护：解析为绝对路径并验证
        self.base_dir = self.base_dir.resolve()
        if ".." in str(self.base_dir):
            raise ValueError("base_dir cannot contain path traversal sequences")

        self.pipeline_name = pipeline_name
        self.checkpoint_file = self.base_dir / f".{pipeline_name}_checkpoint.json"
        self.logger = get_logger()

        # 检查点数据
        self._checkpoints: Dict[str, StepCheckpoint] = {}
        self._current_step: Optional[str] = None

        # 线程锁，确保线程安全
        self._lock = threading.Lock()

        # 加载现有检查点
        self._load()

    def _load(self) -> None:
        """从文件加载检查点"""
        if self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, "r", encoding="utf-8") as f:
                    data = json.load(f)

                if not self._validate_checkpoint_data(data):
                    self.logger.warning(
                        f"Checkpoint file appears corrupted (unexpected interruption). "
                        f"Consider removing {self.checkpoint_file} to start fresh."
                    )
                    self._checkpoints = {}
                    return

                for step_name, step_data in data.get("steps", {}).items():
                    self._checkpoints[step_name] = StepCheckpoint(**step_data)

                self._current_step = data.get("current_step")
                self.logger.info(f"Loaded checkpoint from {self.checkpoint_file}")

                self._check_for_abnormal_interruption()

            except json.JSONDecodeError as e:
                self.logger.warning(
                    f"Checkpoint file is corrupted (JSON parse error: {e}). "
                    f"Please remove {self.checkpoint_file} to start fresh."
                )
                self._checkpoints = {}
            except Exception as e:
                self.logger.warning(f"Failed to load checkpoint: {e}")
                self._checkpoints = {}

    def _validate_checkpoint_data(self, data: Dict[str, Any]) -> bool:
        """验证检查点数据的完整性"""
        if not isinstance(data, dict):
            return False
        if "pipeline_name" not in data:
            return False
        steps = data.get("steps", {})
        if not isinstance(steps, dict):
            return False
        for step_name, step_data in steps.items():
            if not isinstance(step_data, dict):
                return False
            if "status" not in step_data:
                return False
        return True

    def _check_for_abnormal_interruption(self) -> None:
        """检测是否存在异常中断后的损坏状态"""
        for step_name, checkpoint in self._checkpoints.items():
            if checkpoint.status == "running":
                self.logger.warning(
                    f"Checkpoint shows step '{step_name}' was interrupted (status=running). "
                    f"Last run may have been terminated unexpectedly."
                )
                checkpoint.status = "failed"
                checkpoint.error_message = "Interrupted by user or system"
                self.save()

    def save(self) -> None:
        """保存检查点到文件（原子写入）"""
        try:
            data = {
                "pipeline_name": self.pipeline_name,
                "current_step": self._current_step,
                "last_updated": datetime.now().isoformat(),
                "steps": {
                    name: asdict(checkpoint)
                    for name, checkpoint in self._checkpoints.items()
                },
            }

            self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)
            tmp_file = self.checkpoint_file.with_suffix(".tmp")
            with safe_writer(tmp_file, encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            os.chmod(tmp_file, 0o600)
            tmp_file.replace(self.checkpoint_file)

            self.logger.debug(f"Checkpoint saved to {self.checkpoint_file}")

        except Exception as e:
            self.logger.warning(f"Failed to save checkpoint: {e}")
            tmp_file = self.checkpoint_file.with_suffix(".tmp")
            try:
                if tmp_file.exists():
                    tmp_file.unlink()
            except Exception:
                pass

    def start_step(self, step_name: str) -> None:
        """开始执行一个步骤（线程安全）"""
        with self._lock:
            self._current_step = step_name

            if step_name not in self._checkpoints:
                self._checkpoints[step_name] = StepCheckpoint(
                    step_name=step_name,
                    status="running",
                    start_time=datetime.now().isoformat(),
                )
            else:
                self._checkpoints[step_name].status = "running"
                self._checkpoints[step_name].start_time = datetime.now().isoformat()
                self._checkpoints[step_name].error_message = None

            self.save()
            self.logger.info(f"Checkpoint: Started step '{step_name}'")

    def complete_step(
        self, step_name: str, data: Optional[Dict[str, Any]] = None
    ) -> None:
        """完成一个步骤（线程安全）"""
        with self._lock:
            if step_name in self._checkpoints:
                self._checkpoints[step_name].status = "completed"
                self._checkpoints[step_name].end_time = datetime.now().isoformat()

                if data:
                    self._checkpoints[step_name].data.update(data)

                self.save()
                self.logger.info(f"Checkpoint: Completed step '{step_name}'")

    def fail_step(self, step_name: str, error_message: str) -> None:
        """标记步骤失败（线程安全）"""
        with self._lock:
            if step_name in self._checkpoints:
                self._checkpoints[step_name].status = "failed"
                self._checkpoints[step_name].end_time = datetime.now().isoformat()
                self._checkpoints[step_name].error_message = error_message

                self.save()
                self.logger.error(
                    f"Checkpoint: Step '{step_name}' failed - {error_message}"
                )

    def skip_step(self, step_name: str) -> None:
        """跳过步骤（已完成后恢复时）（线程安全）"""
        with self._lock:
            if step_name in self._checkpoints:
                self._checkpoints[step_name].status = "skipped"
                self.save()
                self.logger.info(f"Checkpoint: Skipping completed step '{step_name}'")

    def is_step_completed(self, step_name: str) -> bool:
        """检查步骤是否已完成（线程安全）"""
        with self._lock:
            if step_name in self._checkpoints:
                return self._checkpoints[step_name].status in ("completed", "skipped")
            return False

    def is_step_failed(self, step_name: str) -> bool:
        """检查步骤是否失败（线程安全）"""
        with self._lock:
            if step_name in self._checkpoints:
                return self._checkpoints[step_name].status == "failed"
            return False

    def get_step_data(self, step_name: str) -> Optional[Dict[str, Any]]:
        """获取步骤保存的数据（线程安全）"""
        with self._lock:
            if step_name in self._checkpoints:
                return self._checkpoints[step_name].data
            return None

    def get_step_status(self, step_name: str) -> Optional[str]:
        """获取步骤状态（线程安全）"""
        with self._lock:
            if step_name in self._checkpoints:
                return self._checkpoints[step_name].status
            return None

    def get_all_steps(self) -> List[str]:
        """获取所有已记录的步骤（线程安全）"""
        with self._lock:
            return list(self._checkpoints.keys())

    def reset(self) -> None:
        """重置所有检查点（线程安全）"""
        with self._lock:
            self._checkpoints = {}
            self._current_step = None

            if self.checkpoint_file.exists():
                try:
                    self.checkpoint_file.unlink()
                    self.logger.info(f"Checkpoint file removed: {self.checkpoint_file}")
                except Exception as e:
                    self.logger.warning(f"Failed to remove checkpoint file: {e}")

    def execute_with_checkpoint(
        self,
        step_name: str,
        func: Callable,
        *args: Any,
        skip_if_completed: bool = True,
        **kwargs: Any,
    ) -> Any:
        """
        带检查点保护的执行

        如果步骤已完成且 skip_if_completed=True，则跳过执行
        否则执行函数并在完成后保存检查点
        """
        # 检查是否已完成
        if skip_if_completed and self.is_step_completed(step_name):
            self.skip_step(step_name)
            return self.get_step_data(step_name)

        # 开始步骤
        self.start_step(step_name)

        try:
            # 执行函数
            result = func(*args, **kwargs)

            # 保存结果
            data = {"result": result} if result is not None else {}
            self.complete_step(step_name, data)

            return result

        except Exception as e:
            self.fail_step(step_name, str(e))
            raise


class MCMCTreeCheckpointManager(CheckpointManager):
    """
    MCMCTree 专用检查点管理器

    针对 MCMCTree 的两阶段执行优化
    """

    STEPS = [
        "rate_estimation",  # 速率估计
        "hessian_calculation",  # Hessian 计算
        "mcmc_sampling",  # MCMC 采样
    ]

    def __init__(self, base_dir: Union[str, Path]) -> None:
        super().__init__(base_dir, pipeline_name="mcmctree")

    def get_next_step(self) -> Optional[str]:
        """获取下一个待执行的步骤"""
        for step in self.STEPS:
            if not self.is_step_completed(step):
                return step
        return None

    def get_progress(self) -> Dict[str, Any]:
        """获取执行进度"""
        total = len(self.STEPS)
        completed = sum(1 for step in self.STEPS if self.is_step_completed(step))

        return {
            "total_steps": total,
            "completed_steps": completed,
            "progress_percentage": (completed / total * 100) if total > 0 else 0,
            "current_step": self._current_step,
            "remaining_steps": [s for s in self.STEPS if not self.is_step_completed(s)],
        }


class MethodCheckpoint:
    """单个方法的检查点信息

    ``input_fingerprint`` 记录"这个方法是在哪一次输入指纹下跑完成的"，
    使 ``--resume`` 的跳过判定从"方法名 + 上次是否成功"升级为
    "方法名 + 上次是否成功 + 输入是否还是那份"（审阅项 B-23 的缓存粒度）。
    旧检查点文件里没有这个字段，加载时为 ``None`` —— 按"无法核验"处理。
    """

    def __init__(self, method_name: str) -> None:
        self.method_name = method_name
        self.status: str = "pending"
        self.start_time: Optional[str] = None
        self.end_time: Optional[str] = None
        self.execution_seconds: float = 0.0
        self.result_path: Optional[str] = None
        self.error_message: Optional[str] = None
        self.is_converged: Optional[bool] = None
        self.node_count: int = 0
        self.input_fingerprint: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "method_name": self.method_name,
            "status": self.status,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "execution_seconds": self.execution_seconds,
            "result_path": self.result_path,
            "error_message": self.error_message,
            "is_converged": self.is_converged,
            "node_count": self.node_count,
            "input_fingerprint": self.input_fingerprint,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MethodCheckpoint":
        checkpoint = cls(data["method_name"])
        checkpoint.status = data.get("status", "pending")
        checkpoint.start_time = data.get("start_time")
        checkpoint.end_time = data.get("end_time")
        checkpoint.execution_seconds = data.get("execution_seconds", 0.0)
        checkpoint.result_path = data.get("result_path")
        checkpoint.error_message = data.get("error_message")
        checkpoint.is_converged = data.get("is_converged")
        checkpoint.node_count = data.get("node_count", 0)
        checkpoint.input_fingerprint = data.get("input_fingerprint")
        return checkpoint


class PipelineCheckpointManager:
    """
    管线级检查点管理器

    支持跨方法的断点续跑，记录各方法的执行状态和结果路径。
    当管线中断后重启时，可以跳过已成功完成的方法，
    仅重新执行失败或未执行的方法。

    设计原则：
    - 方法级别的检查点追踪
    - 结果持久化到标准位置
    - 支持从检查点恢复已完成的 DatingResult
    """

    PIPELINE_STEPS = [
        "qc_processing",
        "input_validation",
        "calibration_resolution",
        "method_execution",
        "report_generation",
    ]

    def __init__(self, base_dir: Union[str, Path], max_retries: int = 3) -> None:
        self.base_dir = Path(base_dir) if isinstance(base_dir, str) else base_dir
        self.base_dir = self.base_dir.resolve()

        self.checkpoint_file = self.base_dir / ".phylodater_pipeline_checkpoint.json"
        self.logger = get_logger()
        self.max_retries = max_retries

        self._pipeline_status: str = "pending"
        self._method_checkpoints: Dict[str, MethodCheckpoint] = {}
        self._step_checkpoints: Dict[str, StepCheckpoint] = {}
        self._results_cache: Dict[str, Dict[str, Any]] = {}
        self._failure_counts: Dict[str, int] = {}
        self._input_hashes: Dict[str, str] = {}
        self._run_parameters: Dict[str, Any] = {}

        self._lock = threading.Lock()
        self._load()

    def _load(self) -> None:
        """从文件加载检查点"""
        if not self.checkpoint_file.exists():
            return

        try:
            with open(self.checkpoint_file, "r", encoding="utf-8") as f:
                data = json.load(f)

            self._pipeline_status = data.get("pipeline_status", "pending")

            for method_name, method_data in data.get("methods", {}).items():
                self._method_checkpoints[method_name] = MethodCheckpoint.from_dict(
                    method_data
                )

            for step_name, step_data in data.get("steps", {}).items():
                self._step_checkpoints[step_name] = StepCheckpoint(**step_data)

            self._results_cache = data.get("results_cache", {})
            self._failure_counts = data.get("failure_counts", {})
            self._input_hashes = data.get("input_hashes", {})
            self._run_parameters = data.get("run_parameters", {})

            self.logger.info(f"Loaded pipeline checkpoint from {self.checkpoint_file}")
            self._log_recovery_info()

        except Exception as e:
            self.logger.warning(f"Failed to load pipeline checkpoint: {e}")

    def _log_recovery_info(self) -> None:
        """记录恢复信息"""
        completed_methods = [
            name
            for name, ckpt in self._method_checkpoints.items()
            if ckpt.status == "completed"
        ]
        failed_methods = [
            name
            for name, ckpt in self._method_checkpoints.items()
            if ckpt.status == "failed"
        ]
        pending_methods = [
            name
            for name, ckpt in self._method_checkpoints.items()
            if ckpt.status == "pending"
        ]

        if completed_methods:
            self.logger.info(
                f"Checkpoint: {len(completed_methods)} methods already completed: {completed_methods}"
            )
        if failed_methods:
            self.logger.warning(
                f"Checkpoint: {len(failed_methods)} methods failed: {failed_methods}"
            )
        if pending_methods:
            self.logger.info(
                f"Checkpoint: {len(pending_methods)} methods pending: {pending_methods}"
            )
        if completed_methods and self.get_input_fingerprint() is None:
            # B-23：旧检查点（或从未登记指纹的运行）没有可核验的输入指纹。
            # 复用这些结果在科学上不安全，必须在这里就说清楚。
            self.logger.warning(
                f"Checkpoint has {len(completed_methods)} completed method result(s) but "
                "records no input fingerprint: --resume cannot verify that the tree / "
                "calibration file / configuration are unchanged. The pipeline will "
                "re-run those methods unless record_input_hash(...) is fed a run "
                "fingerprint (see compute_input_hash)."
            )

    def _save_unlocked(self) -> None:
        """保存检查点到文件（调用方必须已持有 self._lock）

        如果磁盘写入失败，回滚内存中的状态变更，
        避免内存状态与磁盘不一致导致断点续跑行为异常。
        """
        try:
            data = {
                "pipeline_status": self._pipeline_status,
                "last_updated": datetime.now().isoformat(),
                "methods": {
                    name: ckpt.to_dict()
                    for name, ckpt in self._method_checkpoints.items()
                },
                "steps": {
                    name: asdict(step) for name, step in self._step_checkpoints.items()
                },
                "results_cache": self._results_cache,
                "failure_counts": self._failure_counts,
                "input_hashes": self._input_hashes,
                "run_parameters": self._run_parameters,
            }

            self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)
            tmp_file = self.checkpoint_file.with_suffix(".tmp")
            with safe_writer(tmp_file, encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            os.chmod(tmp_file, 0o600)
            tmp_file.replace(self.checkpoint_file)

        except Exception as e:
            self.logger.warning(f"Failed to save pipeline checkpoint: {e}")
            # 清理临时文件
            tmp_file = self.checkpoint_file.with_suffix(".tmp")
            try:
                if tmp_file.exists():
                    tmp_file.unlink()
            except Exception:
                pass

    def save(self) -> None:
        """保存检查点到文件（线程安全）"""
        with self._lock:
            self._save_unlocked()

    def start_method(self, method_name: str) -> None:
        """开始执行一个方法"""
        with self._lock:
            if method_name not in self._method_checkpoints:
                self._method_checkpoints[method_name] = MethodCheckpoint(method_name)

            ckpt = self._method_checkpoints[method_name]
            ckpt.status = "running"
            ckpt.start_time = datetime.now().isoformat()
            ckpt.error_message = None

            self._save_unlocked()
            self.logger.info(f"Pipeline checkpoint: Started method '{method_name}'")

    def complete_method(
        self,
        method_name: str,
        result: "DatingResult",
        input_fingerprint: Optional[str] = None,
    ) -> None:
        """完成一个方法。

        ``input_fingerprint`` 是"本次成功结果所对应的输入指纹"（B-23 的缓存粒度）。
        不传时沿用登记在 ``input_hashes`` 里的运行级指纹，因此调用方只需在算出
        指纹的那一处传一次即可。
        """
        with self._lock:
            fingerprint = input_fingerprint or self._current_fingerprint_unlocked()
            if method_name in self._method_checkpoints:
                ckpt = self._method_checkpoints[method_name]
                ckpt.status = "completed"
                ckpt.end_time = datetime.now().isoformat()
                ckpt.execution_seconds = result.execution_seconds
                ckpt.is_converged = result.is_converged
                ckpt.node_count = len(result.node_ages)
                ckpt.input_fingerprint = fingerprint

                self._results_cache[method_name] = result.to_dict()

            self._save_unlocked()
            self.logger.info(f"Pipeline checkpoint: Completed method '{method_name}'")

    def _current_fingerprint_unlocked(self) -> Optional[str]:
        """调用方已持有 ``self._lock`` 时读取运行级输入指纹。"""
        for key in _FINGERPRINT_KEYS:
            value = self._input_hashes.get(key)
            if isinstance(value, str) and value:
                return value
        params_value = self._run_parameters.get("run_fingerprint")
        return params_value if isinstance(params_value, str) and params_value else None

    def fail_method(self, method_name: str, error_message: str) -> None:
        """标记方法失败，记录失败次数"""
        with self._lock:
            if method_name in self._method_checkpoints:
                ckpt = self._method_checkpoints[method_name]
                ckpt.status = "failed"
                ckpt.end_time = datetime.now().isoformat()
                ckpt.error_message = error_message

            self._failure_counts[method_name] = (
                self._failure_counts.get(method_name, 0) + 1
            )

            self._save_unlocked()
            count = self._failure_counts[method_name]
            self.logger.error(
                f"Pipeline checkpoint: Method '{method_name}' failed "
                f"(attempt {count}/{self.max_retries}) - {error_message}"
            )

    def is_method_completed(
        self, method_name: str, input_fingerprint: Optional[str] = None
    ) -> bool:
        """检查方法是否已完成（线程安全）。

        传入 ``input_fingerprint`` 时，"已完成"还要求该方法上次成功时的输入指纹
        与之相同；不同或旧检查点里没有指纹记录 → ``False``（该方法需要重跑）。
        """
        with self._lock:
            ckpt = self._method_checkpoints.get(method_name)
            if ckpt is None or ckpt.status != "completed":
                return False
            if input_fingerprint is None:
                return True
            return bool(
                ckpt.input_fingerprint and ckpt.input_fingerprint == input_fingerprint
            )

    def is_method_failed(self, method_name: str) -> bool:
        """检查方法是否失败（线程安全）"""
        with self._lock:
            return (
                method_name in self._method_checkpoints
                and self._method_checkpoints[method_name].status == "failed"
            )

    def get_cached_result(
        self, method_name: str, input_fingerprint: Optional[str] = None
    ) -> Optional["DatingResult"]:
        """获取缓存的结果用于恢复（线程安全）。

        传入 ``input_fingerprint`` 时，只有缓存结果确实是在**同一份输入**下产出
        的才返回；否则返回 ``None`` 并记 warning（B-23：换树/改校准后不得把上一
        次的年龄当本次结果）。
        """
        with self._lock:
            if method_name not in self._results_cache:
                return None
            cached = self._results_cache[method_name]
            ckpt = self._method_checkpoints.get(method_name)
            stored_fingerprint = ckpt.input_fingerprint if ckpt else None

        if input_fingerprint is not None:
            if not stored_fingerprint:
                self.logger.warning(
                    f"Refusing cached result for '{method_name}': the checkpoint records "
                    "no input fingerprint, so reuse cannot be verified (B-23). "
                    "Re-running this method instead."
                )
                return None
            if stored_fingerprint != input_fingerprint:
                self.logger.warning(
                    f"Refusing cached result for '{method_name}': its input fingerprint "
                    f"({stored_fingerprint[:16]}…) differs from this run's "
                    f"({input_fingerprint[:16]}…). The cached ages belong to another "
                    "input (B-23)."
                )
                return None

        try:
            from ..models import DatingResult

            return DatingResult.from_dict(cached)
        except Exception as e:
            self.logger.warning(
                f"Failed to restore cached result for {method_name}: {e}"
            )
            return None

    def get_methods_to_run(
        self,
        requested_methods: List[str],
        skip_completed: bool = True,
        input_fingerprint: Optional[str] = None,
    ) -> List[str]:
        """获取需要执行的方法列表（跳过已完成的和超过重试上限的）（线程安全）

        Args:
            requested_methods: 请求执行的方法列表
            skip_completed: 是否跳过已完成的方法（默认 True，--resume 时为 True）
            input_fingerprint: 本次运行的输入指纹。提供时，"已完成"还必须满足
                "上次成功时的输入指纹与之相同"，否则不得跳过（B-23：仅凭方法名 +
                上次成功就把换过输入的运行直接复用，会把上一次的年龄当本次结果）。
                不提供时保持旧行为——是否复用由调用方自行核验。
        """
        with self._lock:
            methods_to_run = []
            for method_name in requested_methods:
                ckpt = self._method_checkpoints.get(method_name)
                if (
                    skip_completed
                    and ckpt is not None
                    and ckpt.status == "completed"
                    and self._fingerprint_allows_reuse(ckpt, input_fingerprint)
                ):
                    self.logger.info(f"Skipping completed method: {method_name}")
                elif ckpt is not None and ckpt.status == "failed":
                    failure_count = self._failure_counts.get(method_name, 0)
                    if failure_count >= self.max_retries:
                        self.logger.warning(
                            f"Skipping method '{method_name}': exceeded max retries "
                            f"({failure_count}/{self.max_retries})"
                        )
                    else:
                        self.logger.warning(f"Will retry failed method: {method_name}")
                        methods_to_run.append(method_name)
                else:
                    if (
                        skip_completed
                        and ckpt is not None
                        and ckpt.status == "completed"
                    ):
                        self.logger.warning(
                            f"Will re-run completed method '{method_name}': 输入指纹与"
                            "检查点记录不一致或检查点无指纹可核验（B-23）"
                        )
                    methods_to_run.append(method_name)

            return methods_to_run

    @staticmethod
    def _fingerprint_allows_reuse(
        ckpt: "MethodCheckpoint", input_fingerprint: Optional[str]
    ) -> bool:
        """（持锁调用）方法级复用判定：没传指纹时一律允许复用（旧行为）。"""
        if input_fingerprint is None:
            return True
        return bool(
            ckpt.input_fingerprint and ckpt.input_fingerprint == input_fingerprint
        )

    def start_pipeline_step(self, step_name: str) -> None:
        """开始一个管线步骤"""
        with self._lock:
            if step_name not in self._step_checkpoints:
                self._step_checkpoints[step_name] = StepCheckpoint(
                    step_name=step_name,
                    status="running",
                    start_time=datetime.now().isoformat(),
                )
            else:
                self._step_checkpoints[step_name].status = "running"
                self._step_checkpoints[step_name].start_time = (
                    datetime.now().isoformat()
                )

            self._save_unlocked()

    def complete_pipeline_step(self, step_name: str) -> None:
        """完成一个管线步骤"""
        with self._lock:
            if step_name in self._step_checkpoints:
                self._step_checkpoints[step_name].status = "completed"
                self._step_checkpoints[step_name].end_time = datetime.now().isoformat()

            self._save_unlocked()

    def is_step_completed(self, step_name: str) -> bool:
        """检查步骤是否已完成"""
        return step_name in self._step_checkpoints and self._step_checkpoints[
            step_name
        ].status in ("completed", "skipped")

    def get_pipeline_status(self) -> str:
        """获取管线状态"""
        return self._pipeline_status

    def set_pipeline_status(self, status: str) -> None:
        """设置管线状态"""
        self._pipeline_status = status
        self.save()

    def record_input_hash(self, key: Union[str, Path], hash_value: str) -> None:
        """记录一份输入（或运行级指纹）的哈希值（线程安全）。

        ``key`` 既接受真实文件路径（``Path``），也接受**逻辑键**——例如
        ``RUN_FINGERPRINT_KEY``（"phylodater:run-fingerprint"），后者用来存整次
        运行的聚合输入指纹（B-23）。内部一律以 ``str(key)`` 存储。

        Raises:
            TypeError: ``key`` 不是字符串/Path，或 ``hash_value`` 不是字符串。
            ValueError: ``hash_value`` 为空——空指纹会让 ``input_matches`` 误判为
                "已核验一致"，宁可在登记时就失败。
        """
        normalized_key = self._normalize_hash_key(key)
        if not isinstance(hash_value, str):
            raise TypeError(
                f"hash_value must be a str digest, got {type(hash_value).__name__}"
            )
        if not hash_value.strip():
            raise ValueError("hash_value must be a non-empty digest string")
        with self._lock:
            self._input_hashes[normalized_key] = hash_value
            self._save_unlocked()

    def get_input_hash(self, key: Union[str, Path]) -> Optional[str]:
        """读取此前登记的某个输入哈希（B-23 缺的那半个"读"端）。

        从未登记过该键（含因输入变化而重置）时返回 ``None``——调用方必须把
        ``None`` 当作"无从确认是否一致"，按"已变更"处理，而不是当作"一致"。
        """
        try:
            normalized_key = self._normalize_hash_key(key)
        except TypeError:
            return None
        with self._lock:
            value = self._input_hashes.get(normalized_key)
        return value if isinstance(value, str) else None

    def record_run_parameters(self, params: Dict[str, Any]) -> None:
        """记录运行参数（线程安全）。

        Raises:
            TypeError: ``params`` 不是映射。
        """
        if not isinstance(params, Mapping):
            raise TypeError(f"params must be a mapping, got {type(params).__name__}")
        with self._lock:
            self._run_parameters.update(dict(params))
            self._save_unlocked()

    @staticmethod
    def _normalize_hash_key(key: Union[str, Path]) -> str:
        """把逻辑键 / 文件路径统一成字符串键。"""
        if isinstance(key, Path):
            return str(key)
        if isinstance(key, str):
            if not key.strip():
                raise ValueError("input hash key must be a non-empty string")
            return key
        raise TypeError(
            f"input hash key must be a str or Path, got {type(key).__name__}"
        )

    @staticmethod
    def compute_input_hash(inputs: Mapping[Any, Any]) -> str:
        """由规范化输入映射算出稳定的输入指纹字符串（B-23）。

        与 ``core/pipeline.py`` 的本地兜底实现保持**同一口径**
        （``sha256:`` + 对按键排序的 JSON 规范化串取哈希），因此两侧写出的指纹
        可以直接互比，不会因为换了实现就把一次本可复用的运行判成"输入已变"。

        注意：本函数只做**确定性规范化**，不读文件；调用方负责把树内容、校准
        文件内容、生效配置、外部软件版本等已经哈希/字符串化的量传进来。
        """
        if not isinstance(inputs, Mapping):
            raise TypeError(
                f"compute_input_hash expects a mapping, got {type(inputs).__name__}"
            )
        import hashlib

        normalized = {str(key): value for key, value in inputs.items()}
        canonical = json.dumps(normalized, sort_keys=True, default=repr)
        return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def get_input_fingerprint(self) -> Optional[str]:
        """检查点里登记的"上一次运行的输入指纹"，没有则 ``None``。

        取用顺位：逻辑键 ``RUN_FINGERPRINT_KEY`` 等专用键 → ``run_parameters``
        里的 ``run_fingerprint`` → 由逐文件哈希现算的聚合指纹。
        """
        with self._lock:
            for key in _FINGERPRINT_KEYS:
                value = self._input_hashes.get(key)
                if isinstance(value, str) and value:
                    return value
            params_value = self._run_parameters.get("run_fingerprint")
            if isinstance(params_value, str) and params_value:
                return params_value
            individuals = {
                key: value
                for key, value in self._input_hashes.items()
                if key not in _FINGERPRINT_KEYS
            }
        if individuals:
            return self.compute_input_hash(individuals)
        return None

    def input_matches(self, current_hash: str) -> bool:
        """本次输入指纹是否与检查点里记录的一致（B-23 的判定入口）。

        语义（严格 fail-closed）：
          * 只有"检查点里存过指纹 **且** 与传入指纹逐字相等"才返回 ``True``；
          * 没有存过任何指纹（旧版本检查点、或从没登记过）→ ``False``：
            无从判断输入是否变过，绝不能把上一次的年龄当本次结果复用；
          * 传入指纹为空/非字符串 → 同样 ``False``。
        """
        if not isinstance(current_hash, str) or not current_hash.strip():
            self.logger.info(
                "input_matches(): 本次输入指纹为空，按'无法核验'处理 → 不复用缓存结果"
            )
            return False
        stored = self.get_input_fingerprint()
        if stored is None:
            self.logger.info(
                "input_matches(): 检查点中没有输入指纹记录，按'输入可能已变更'处理 "
                "→ 不复用缓存结果（请用 record_input_hash(...) 登记运行指纹）"
            )
            return False
        matches = stored == current_hash
        if not matches:
            self.logger.warning(
                f"输入指纹与检查点不一致（本次 {current_hash[:16]}… ≠ "
                f"检查点 {stored[:16]}…），已缓存结果对应另一次输入"
            )
        return matches

    def check_input_changed(self, key: Union[str, Path], hash_value: str) -> bool:
        """检查某个输入是否已变更（没有记录时同样返回 ``True`` = 已变更）。"""
        try:
            old_hash = self._input_hashes.get(self._normalize_hash_key(key))
        except (TypeError, ValueError):
            return True
        if old_hash is None:
            return True
        return old_hash != hash_value

    def reset(self) -> None:
        """重置所有检查点并删除检查点文件。

        B-23：``_input_hashes`` / ``_run_parameters`` 必须一起清空——否则上一轮
        留下的输入指纹会"复活"，让下一次运行在**真正换过输入**的情况下被判为
        "指纹一致"，从而复用根本不该复用的陈旧结果。
        """
        with self._lock:
            self._pipeline_status = "pending"
            self._method_checkpoints = {}
            self._step_checkpoints = {}
            self._results_cache = {}
            self._failure_counts = {}
            self._input_hashes = {}
            self._run_parameters = {}

            if self.checkpoint_file.exists():
                try:
                    self.checkpoint_file.unlink()
                    self.logger.info("Pipeline checkpoint file removed")
                except Exception as e:
                    self.logger.warning(f"Failed to remove checkpoint file: {e}")

    def get_progress_summary(self) -> Dict[str, Any]:
        """获取进度摘要（线程安全）"""
        with self._lock:
            total_methods = len(self._method_checkpoints)
            completed_methods = sum(
                1
                for ckpt in self._method_checkpoints.values()
                if ckpt.status == "completed"
            )
            failed_methods = sum(
                1
                for ckpt in self._method_checkpoints.values()
                if ckpt.status == "failed"
            )
            fingerprint = self._current_fingerprint_unlocked()

            return {
                "pipeline_status": self._pipeline_status,
                "total_methods": total_methods,
                "completed_methods": completed_methods,
                "failed_methods": failed_methods,
                "pending_methods": total_methods - completed_methods - failed_methods,
                "max_retries": self.max_retries,
                "failure_counts": dict(self._failure_counts),
                # B-23：让"这次续跑到底有没有输入可核验"成为一个可查的事实
                "input_fingerprint": fingerprint,
                "resume_reuse_verifiable": fingerprint is not None,
                "methods_without_fingerprint": [
                    name
                    for name, ckpt in self._method_checkpoints.items()
                    if ckpt.status == "completed" and not ckpt.input_fingerprint
                ],
                "methods": {
                    name: ckpt.to_dict()
                    for name, ckpt in self._method_checkpoints.items()
                },
            }
