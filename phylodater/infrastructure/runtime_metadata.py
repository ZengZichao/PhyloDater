"""
Runtime Metadata - 机器可读运行元数据

生成标准化的运行时元数据 JSON 文件，用于：
- 计算可重复性记录
- 资源使用审计
- 执行过程追溯
"""

import json
import os
import platform
import threading
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

from .safe_io import safe_writer


class RunStatus(Enum):
    """运行状态"""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    PARTIAL = "partial"  # 部分方法成功
    INTERRUPTED = "interrupted"  # 用户中断


class StepStatus(Enum):
    """步骤状态"""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class MethodMetadata:
    """单个方法的元数据"""

    status: str
    execution_seconds: float = 0.0
    is_converged: Optional[bool] = None
    error_message: Optional[str] = None
    steps: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    software_versions: Dict[str, str] = field(default_factory=dict)
    node_count: int = 0
    warnings: List[str] = field(default_factory=list)


@dataclass
class RuntimeMetadata:
    """
    运行时元数据容器

    生成符合以下标准的元数据文件：
    - pipeline: 管线名称
    - version: 管线版本
    - run_id: 唯一运行标识符
    - start_time/end_time: ISO 8601 格式时间戳
    - duration_seconds: 总运行时长
    - status: 运行状态
    - seed: 随机种子（用于可重现性）
    - methods: 各方法的执行元数据
    - resources: 系统资源信息
    """

    pipeline: str = "phylodater"
    version: str = "0.1.0"
    run_id: str = ""
    start_time: str = ""
    end_time: str = ""
    duration_seconds: float = 0.0
    status: str = RunStatus.PENDING.value
    seed: Optional[int] = None
    methods: Dict[str, MethodMetadata] = field(default_factory=dict)
    resources: Dict[str, Any] = field(default_factory=dict)
    input_files: Dict[str, str] = field(default_factory=dict)
    output_dir: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典用于 JSON 序列化"""
        # 必须显式标成 Dict[str, Any]：字面量里既有 str/float 也有嵌套 dict，
        # 不标注时 mypy 推不出 ``result["methods"]`` 还能按下标写。
        result: Dict[str, Any] = {
            "pipeline": self.pipeline,
            "version": self.version,
            "run_id": self.run_id,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "duration_seconds": self.duration_seconds,
            "status": self.status,
            "seed": self.seed,
            "methods": {},
            "resources": self.resources,
            "input_files": self.input_files,
            "output_dir": self.output_dir,
        }

        for name, method_meta in self.methods.items():
            result["methods"][name] = asdict(method_meta)

        return result

    def to_json(self, output_path: Path) -> None:
        """保存为 JSON 文件"""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with safe_writer(output_path, encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RuntimeMetadata":
        """从字典创建"""
        metadata = cls(
            pipeline=data.get("pipeline", "phylodater"),
            version=data.get("version", "0.1.0"),
            run_id=data.get("run_id", ""),
            start_time=data.get("start_time", ""),
            end_time=data.get("end_time", ""),
            duration_seconds=data.get("duration_seconds", 0.0),
            status=data.get("status", RunStatus.PENDING.value),
            seed=data.get("seed"),
            resources=data.get("resources", {}),
            input_files=data.get("input_files", {}),
            output_dir=data.get("output_dir", ""),
        )

        for name, method_data in data.get("methods", {}).items():
            metadata.methods[name] = MethodMetadata(**method_data)

        return metadata

    @classmethod
    def from_json(cls, input_path: Path) -> "RuntimeMetadata":
        """从 JSON 文件加载"""
        with open(input_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return cls.from_dict(data)


class RuntimeMetadataManager:
    """
    运行时元数据管理器

    在管线执行过程中收集和记录元数据
    """

    def __init__(
        self, output_dir: Path, pipeline: str = "phylodater", version: str = "0.1.0"
    ) -> None:
        self.output_dir = output_dir
        self.metadata = RuntimeMetadata(
            pipeline=pipeline,
            version=version,
            run_id=self._generate_run_id(),
            start_time=self._get_utc_timestamp(),
            output_dir=str(output_dir),
        )
        self._start_wall_time = time.time()
        self._method_start_times: Dict[str, float] = {}
        self._lock = threading.Lock()

    def _generate_run_id(self) -> str:
        """生成唯一运行标识符（使用独立随机源，不污染全局状态）"""
        import secrets

        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        short_id = secrets.token_hex(3)  # 6 个十六进制字符
        return f"{timestamp}_{short_id}"

    def _get_utc_timestamp(self) -> str:
        """获取当前 UTC 时间戳（ISO 8601 格式）"""
        return datetime.now(timezone.utc).isoformat()

    def set_seed(self, seed: Optional[int]) -> None:
        """记录本次运行的随机种子；未固定种子时写 None。

        ``RuntimeMetadata.seed`` 本身就是 ``Optional[int]``，而 pipeline 的
        ``config.seed`` 也允许 None（"不固定种子"是合法运行形态），所以这里
        必须接 Optional[int]，否则调用方只能自己判空。
        """
        self.metadata.seed = seed

    def set_input_file(self, name: str, path: str) -> None:
        """记录输入文件"""
        self.metadata.input_files[name] = path

    def collect_system_resources(self) -> None:
        """收集系统资源信息。

        ``psutil`` 在 ``pyproject.toml`` 里是**可选**依赖（仅用于内存审计），
        此前这里是函数体顶部的硬 ``import psutil``：模块没装时
        ``ModuleNotFoundError`` 会直接穿过 ``try`` 块抛出，把整条流水线在
        ``_init_metadata_manager()`` 处炸掉（pipeline 集成测试与批量测试都
        因此失败）。现在按可选依赖处理：缺 psutil 时降级为"只记录平台信息"，
        并在元数据里显式写明缺失原因，而不是谎报资源数据。
        """
        resources: Dict[str, Any] = {
            "platform": platform.platform(),
            "python_version": platform.python_version(),
            "cpu_count": os.cpu_count() or 1,
            "hostname": platform.node(),
        }

        try:
            import psutil
        except ImportError as e:
            resources["cpu_count_physical"] = None
            resources["memory_total_mb"] = None
            resources["memory_available_mb"] = None
            resources["resource_note"] = (
                f"psutil unavailable ({e.__class__.__name__}: {e}); "
                "memory metrics not recorded and cpu_count is os.cpu_count()"
            )
        else:
            try:
                memory = psutil.virtual_memory()
                resources.update(
                    {
                        "cpu_count": psutil.cpu_count(logical=True)
                        or resources["cpu_count"],
                        "cpu_count_physical": psutil.cpu_count(logical=False),
                        "memory_total_mb": memory.total / (1024 * 1024),
                        "memory_available_mb": memory.available / (1024 * 1024),
                    }
                )
            except (OSError, NotImplementedError, ValueError) as e:
                # 受限容器/沙箱里 psutil 的系统调用会失败（NotImplementedError
                # 是其文档化的行为之一），此时只降级资源审计，不影响定年结果。
                resources["cpu_count_physical"] = None
                resources["memory_total_mb"] = None
                resources["memory_available_mb"] = None
                resources["resource_note"] = (
                    f"psutil present but system memory query failed "
                    f"({e.__class__.__name__}: {e}); memory metrics not recorded"
                )

        self.metadata.resources = resources

    def start_method(self, method_name: str) -> None:
        """记录方法开始执行（线程安全）"""
        with self._lock:
            self._method_start_times[method_name] = time.time()
            if method_name not in self.metadata.methods:
                self.metadata.methods[method_name] = MethodMetadata(
                    status=StepStatus.RUNNING.value
                )
            else:
                self.metadata.methods[method_name].status = StepStatus.RUNNING.value

    def complete_method(
        self,
        method_name: str,
        execution_seconds: float,
        is_converged: Optional[bool] = None,
        node_count: int = 0,
        warnings: Optional[List[str]] = None,
    ) -> None:
        """记录方法完成（线程安全）"""
        with self._lock:
            if method_name in self.metadata.methods:
                method_meta = self.metadata.methods[method_name]
                method_meta.status = StepStatus.COMPLETED.value
                method_meta.execution_seconds = execution_seconds
                method_meta.is_converged = is_converged
                method_meta.node_count = node_count
                if warnings:
                    method_meta.warnings = warnings

    def fail_method(self, method_name: str, error_message: str) -> None:
        """记录方法失败（线程安全）"""
        with self._lock:
            if method_name in self.metadata.methods:
                method_meta = self.metadata.methods[method_name]
                method_meta.status = StepStatus.FAILED.value
                method_meta.error_message = error_message

    def set_method_version(
        self, method_name: str, tool_name: str, version: str
    ) -> None:
        """记录方法使用的软件版本"""
        if method_name not in self.metadata.methods:
            self.metadata.methods[method_name] = MethodMetadata(
                status=StepStatus.PENDING.value
            )

        method_meta = self.metadata.methods[method_name]
        if not method_meta.software_versions:
            method_meta.software_versions = {}
        method_meta.software_versions[tool_name] = version

    def add_step(
        self,
        method_name: str,
        step_name: str,
        status: str,
        duration_seconds: Optional[float] = None,
        details: Optional[Dict[str, Any]] = None,
    ) -> None:
        """添加方法内的步骤信息"""
        if method_name not in self.metadata.methods:
            self.metadata.methods[method_name] = MethodMetadata(
                status=StepStatus.RUNNING.value
            )

        step_info: Dict[str, Any] = {"status": status}
        if duration_seconds is not None:
            step_info["duration_seconds"] = duration_seconds
        if details:
            step_info.update(details)

        self.metadata.methods[method_name].steps[step_name] = step_info

    def finalize(self, status: Optional[RunStatus] = None) -> RuntimeMetadata:
        """完成元数据收集

        Args:
            status: 显式指定的运行状态。如果为 None，则根据方法完成情况自动推断。
        """
        self.metadata.end_time = self._get_utc_timestamp()
        self.metadata.duration_seconds = time.time() - self._start_wall_time

        if status is not None:
            self.metadata.status = status.value
        else:
            # 自动推断状态
            completed_methods = sum(
                1
                for m in self.metadata.methods.values()
                if m.status == StepStatus.COMPLETED.value
            )
            total_methods = len(self.metadata.methods)

            if completed_methods == 0:
                self.metadata.status = RunStatus.FAILED.value
            elif completed_methods < total_methods:
                self.metadata.status = RunStatus.PARTIAL.value
            else:
                self.metadata.status = RunStatus.COMPLETED.value

        return self.metadata

    def save(self, filename: str = "runtime_metadata.json") -> Path:
        """保存元数据到文件"""
        output_path = self.output_dir / filename
        self.metadata.to_json(output_path)
        return output_path


def load_metadata(
    output_dir: Path, filename: str = "runtime_metadata.json"
) -> RuntimeMetadata:
    """加载已保存的运行时元数据"""
    return RuntimeMetadata.from_json(output_dir / filename)
