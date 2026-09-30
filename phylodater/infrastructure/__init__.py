"""
PhyloDater 基础设施层

提供日志、进程管理、配置管理等基础服务
"""

from .alignment_metadata import AlignmentMetadata, AlignmentMetadataExtractor
from .auto_naming import AutoNamingManager, LSD2AutoNaming, MCMCTreeAutoNaming
from .checkpoint import (
    CheckpointManager,
    CheckpointStatus,
    MCMCTreeCheckpointManager,
    StepCheckpoint,
)
from .configuration import (
    Configuration,
    LSD2Config,
    MCMCTreeConfig,
    MDCatConfig,
    PATHd8Config,
    R8sConfig,
    SoftwarePaths,
    ToolConfig,
    TreePLConfig,
    WLogDateConfig,
)
from .container import (
    ContainerEnvironment,
    ContainerType,
    create_container_compatible_workdir,
    detect_container_type,
    get_container_environment_info,
    get_writable_directory,
    is_containerized,
)
from .logging import LogLevel, PhyloDaterLogger, get_logger, setup_logging
from .plugins import (
    Plugin,
    PluginLoader,
    PluginMetadata,
    PluginRegistry,
    PluginType,
    QCPlugin,
    QCPluginChain,
    get_plugin_loader,
    get_plugin_registry,
    load_all_plugins,
    register_qc_plugin,
)
from .process_runner import ProcessResult, ProcessRunner
from .runtime_metadata import (
    MethodMetadata,
    RunStatus,
    RuntimeMetadata,
    RuntimeMetadataManager,
    StepStatus,
    load_metadata,
)
from .signal_handler import (
    GracefulShutdownHandler,
    ShutdownContext,
    get_shutdown_handler,
)

__all__ = [
    # 日志
    "PhyloDaterLogger",
    "get_logger",
    "LogLevel",
    "setup_logging",
    # 进程管理
    "ProcessRunner",
    "ProcessResult",
    # 配置
    "Configuration",
    "ToolConfig",
    "SoftwarePaths",
    "MCMCTreeConfig",
    "R8sConfig",
    "TreePLConfig",
    "PATHd8Config",
    "LSD2Config",
    "WLogDateConfig",
    "MDCatConfig",
    # 比对元数据
    "AlignmentMetadataExtractor",
    "AlignmentMetadata",
    # 检查点
    "CheckpointManager",
    "MCMCTreeCheckpointManager",
    "CheckpointStatus",
    "StepCheckpoint",
    # 自动命名
    "AutoNamingManager",
    "LSD2AutoNaming",
    "MCMCTreeAutoNaming",
    # 运行时元数据
    "RuntimeMetadata",
    "RuntimeMetadataManager",
    "RunStatus",
    "StepStatus",
    "MethodMetadata",
    "load_metadata",
    # 容器化兼容
    "ContainerType",
    "ContainerEnvironment",
    "detect_container_type",
    "get_writable_directory",
    "get_container_environment_info",
    "is_containerized",
    "create_container_compatible_workdir",
    # 插件系统
    "PluginType",
    "PluginMetadata",
    "Plugin",
    "QCPlugin",
    "PluginRegistry",
    "PluginLoader",
    "QCPluginChain",
    "get_plugin_registry",
    "get_plugin_loader",
    "load_all_plugins",
    "register_qc_plugin",
    # 信号处理
    "GracefulShutdownHandler",
    "ShutdownContext",
    "get_shutdown_handler",
]
