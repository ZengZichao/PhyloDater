"""
Plugin System - 插件加载机制

支持以插件形式动态加载以下组件：
- 质控步骤 (QC plugins)
- 自定义定年方法 (DatingMethod plugins)
- 自定义可视化后端 (PlotterBackend plugins)

插件发现机制：
1. 从配置指定目录加载
2. 从环境变量 PHYLODATER_PLUGINS 加载
3. 从 ~/.phylodater/plugins 加载
"""

import importlib
import importlib.util
import os
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Type

from .logging import get_logger

logger = get_logger()


class PluginType(Enum):
    """插件类型"""

    QC = "qc"  # 质控步骤
    DATING_METHOD = "dating"  # 定年方法
    VISUALIZATION = "viz"  # 可视化后端
    CONSTRAINT = "constraint"  # 自定义约束类型


@dataclass
class PluginMetadata:
    """插件元数据"""

    name: str
    version: str
    description: str
    plugin_type: PluginType
    author: Optional[str] = None
    license: Optional[str] = None
    dependencies: List[str] = field(default_factory=list)
    config_schema: Optional[Dict[str, Any]] = None


class Plugin(ABC):
    """
    插件抽象基类

    所有插件必须继承此类并实现相应接口
    """

    @abstractmethod
    def get_metadata(self) -> PluginMetadata:
        """返回插件元数据"""
        pass

    @abstractmethod
    def initialize(self, config: Dict[str, Any]) -> None:
        """初始化插件

        Args:
            config: 插件配置
        """
        pass

    def validate_config(self, config: Dict[str, Any]) -> bool:
        """验证配置是否有效

        Args:
            config: 待验证的配置

        Returns:
            配置是否有效
        """
        return True


class QCPlugin(Plugin):
    """
    质控插件基类

    用于在定年前/后对数据进行处理
    """

    @abstractmethod
    def process_tree(self, tree: Any, context: Dict[str, Any]) -> Any:
        """处理树对象

        Args:
            tree: 输入的系统发育树
            context: 上下文信息，包含 alignment、calibrations 等

        Returns:
            处理后的树对象
        """
        pass

    def process_alignment(self, alignment: Any, context: Dict[str, Any]) -> Any:
        """处理比对文件（可选实现）

        Args:
            alignment: 输入的比对对象
            context: 上下文信息

        Returns:
            处理后的比对对象
        """
        return alignment


class DatingMethodPlugin(Plugin):
    """
    定年方法插件基类

    用于添加自定义定年方法
    """

    @abstractmethod
    def create_method(self, config: Any, output_dir: Path) -> Any:
        """创建定年方法实例

        Args:
            config: 工具配置
            output_dir: 输出目录

        Returns:
            DatingMethod 实例
        """
        pass


class VisualizationPlugin(Plugin):
    """
    可视化插件基类

    用于添加自定义可视化后端
    """

    @abstractmethod
    def create_backend(self) -> Any:
        """创建可视化后端实例

        Returns:
            PlotterBackend 实例
        """
        pass


@dataclass
class LoadedPlugin:
    """已加载的插件"""

    metadata: PluginMetadata
    instance: Plugin
    module_path: str
    enabled: bool = True


class PluginRegistry:
    """
    插件注册表

    管理所有已加载的插件
    """

    def __init__(self) -> None:
        self._plugins: Dict[str, LoadedPlugin] = {}
        self._qc_hooks: List[QCPlugin] = []
        self._method_factories: Dict[str, Type[Plugin]] = {}
        self._backend_factories: Dict[str, Type[Plugin]] = {}

    def register_plugin(self, name: str, plugin: LoadedPlugin) -> None:
        """注册插件"""
        self._plugins[name] = plugin
        logger.debug(f"Registered plugin: {name} ({plugin.metadata.plugin_type.value})")

    def get_plugin(self, name: str) -> Optional[LoadedPlugin]:
        """获取插件"""
        return self._plugins.get(name)

    def list_plugins(
        self, plugin_type: Optional[PluginType] = None
    ) -> List[PluginMetadata]:
        """列出插件"""
        plugins = []
        for plugin in self._plugins.values():
            if plugin_type is None or plugin.metadata.plugin_type == plugin_type:
                plugins.append(plugin.metadata)
        return plugins

    def is_enabled(self, name: str) -> bool:
        """检查插件是否启用"""
        plugin = self._plugins.get(name)
        return plugin.enabled if plugin else False

    def enable(self, name: str) -> bool:
        """启用插件"""
        plugin = self._plugins.get(name)
        if plugin:
            plugin.enabled = True
            return True
        return False

    def disable(self, name: str) -> bool:
        """禁用插件"""
        plugin = self._plugins.get(name)
        if plugin:
            plugin.enabled = False
            return True
        return False


class PluginLoader:
    """
    插件加载器

    从多个来源发现和加载插件
    """

    def __init__(self, registry: PluginRegistry) -> None:
        self.registry = registry
        self._discovered_paths: Set[Path] = set()

    def discover_plugins(self, search_paths: Optional[List[Path]] = None) -> List[Path]:
        """
        发现插件

        Args:
            search_paths: 搜索路径列表，None 则使用默认路径

        Returns:
            发现的插件文件路径列表
        """
        if search_paths is None:
            search_paths = self._get_default_search_paths()

        discovered: List[Path] = []

        for path in search_paths:
            if not path.exists():
                continue

            if path.is_file() and path.suffix == ".py":
                discovered.append(path)
            elif path.is_dir():
                for plugin_file in path.glob("*.py"):
                    if not plugin_file.name.startswith("_"):
                        discovered.append(plugin_file)

        self._discovered_paths.update(discovered)
        return list(self._discovered_paths)

    def _get_default_search_paths(self) -> List[Path]:
        """获取默认搜索路径"""
        paths: List[Path] = []

        env_plugins = os.environ.get("PHYLODATER_PLUGINS")
        if env_plugins:
            for p in env_plugins.split(os.pathsep):
                paths.append(Path(p))

        user_plugins = Path.home() / ".phylodater" / "plugins"
        if user_plugins.exists():
            paths.append(user_plugins)

        current_dir = Path.cwd() / "plugins"
        if current_dir.exists():
            paths.append(current_dir)

        return paths

    def load_plugin_from_file(
        self, plugin_path: Path, allowed_dirs: Optional[List[Path]] = None
    ) -> Optional[LoadedPlugin]:
        """
        从文件加载插件

        安全措施：
        - 仅从明确允许的目录加载（不跟随符号链接到敏感目录）
        - 插件文件必须在已知搜索路径内
        - 记录所有插件加载行为用于审计

        Args:
            plugin_path: 插件文件路径

        Returns:
            LoadedPlugin 或 None
        """
        # 安全检查：验证插件路径在已知搜索路径内
        resolved_path = plugin_path.resolve()
        allowed = (
            allowed_dirs
            if allowed_dirs is not None
            else self._get_default_search_paths()
        )
        allowed_dirs_resolved = [p.resolve() for p in allowed if p.exists()]

        # 使用真实祖先判定：allowed_dir 必须是 resolved_path 的某个父目录，
        # 而非前缀字符串匹配（前缀匹配会被 "../" 或同级目录绕过，例如
        # "/plugins_evil" 会以 "/plugins" 开头而误判为允许）。
        in_allowed_dir = any(
            allowed_dir in resolved_path.parents
            for allowed_dir in allowed_dirs_resolved
        )

        if not in_allowed_dir:
            logger.warning(
                f"Rejected plugin from untrusted path: {plugin_path}. "
                f"Plugins must be in: {allowed}"
            )
            return None

        try:
            spec = importlib.util.spec_from_file_location(
                f"phylodater_plugin_{plugin_path.stem}", plugin_path
            )
            if spec is None or spec.loader is None:
                return None

            module = importlib.util.module_from_spec(spec)
            # 使用唯一的模块名避免 sys.modules 冲突
            module_name = f"phylodater_plugin_{resolved_path.stem}_{id(resolved_path)}"
            sys.modules[module_name] = module
            spec.loader.exec_module(module)

            if hasattr(module, "PLUGIN_METADATA"):
                metadata_dict = module.PLUGIN_METADATA
                # 验证必要的元数据字段
                if not isinstance(metadata_dict, dict) or "name" not in metadata_dict:
                    logger.warning(f"Plugin {plugin_path}: invalid PLUGIN_METADATA")
                    return None

                metadata = PluginMetadata(
                    name=metadata_dict["name"],
                    version=metadata_dict.get("version", "0.1.0"),
                    description=metadata_dict.get("description", ""),
                    plugin_type=PluginType(metadata_dict["plugin_type"]),
                    author=metadata_dict.get("author"),
                    license=metadata_dict.get("license"),
                    dependencies=metadata_dict.get("dependencies", []),
                    config_schema=metadata_dict.get("config_schema"),
                )

                if hasattr(module, "PluginClass"):
                    plugin_instance = module.PluginClass()
                    # 仅传递元数据中的 config 字典，不传递任意数据
                    safe_config = metadata_dict.get("config", {})
                    if not isinstance(safe_config, dict):
                        safe_config = {}
                    plugin_instance.initialize(safe_config)

                    logger.info(
                        f"Loaded plugin: {metadata.name} v{metadata.version} from {plugin_path}"
                    )

                    return LoadedPlugin(
                        metadata=metadata,
                        instance=plugin_instance,
                        module_path=str(plugin_path),
                    )

        except Exception as e:
            logger.warning(f"Failed to load plugin from {plugin_path}: {e}")

        return None

    def load_all(self, search_paths: Optional[List[Path]] = None) -> int:
        """
        加载所有发现的插件

        Args:
            search_paths: 搜索路径

        Returns:
            成功加载的插件数量
        """
        plugin_paths = self.discover_plugins(search_paths)
        loaded_count = 0

        for path in plugin_paths:
            plugin = self.load_plugin_from_file(path, allowed_dirs=search_paths)
            if plugin:
                self.registry.register_plugin(plugin.metadata.name, plugin)
                loaded_count += 1

        logger.info(
            f"Loaded {loaded_count} plugins from {len(plugin_paths)} discovered paths"
        )
        return loaded_count


class QCPluginChain:
    """
    质控插件链

    按顺序执行多个质控插件
    """

    def __init__(self) -> None:
        self._plugins: List[QCPlugin] = []

    def add(self, plugin: QCPlugin, position: Optional[int] = None) -> "QCPluginChain":
        """
        添加质控插件

        Args:
            plugin: QCPlugin 实例
            position: 插入位置，None 则追加到末尾

        Returns:
            self (用于链式调用)
        """
        if position is None:
            self._plugins.append(plugin)
        else:
            self._plugins.insert(position, plugin)
        return self

    def process_tree(self, tree: Any, context: Dict[str, Any]) -> Any:
        """
        按顺序处理树

        Args:
            tree: 输入树
            context: 上下文

        Returns:
            处理后的树
        """
        result = tree
        for plugin in self._plugins:
            logger.debug(f"Running QC plugin: {plugin.get_metadata().name}")
            result = plugin.process_tree(result, context)
        return result

    def process_alignment(self, alignment: Any, context: Dict[str, Any]) -> Any:
        """按顺序处理比对"""
        result = alignment
        for plugin in self._plugins:
            if hasattr(plugin, "process_alignment"):
                result = plugin.process_alignment(result, context)
        return result

    def get_plugins(self) -> List[QCPlugin]:
        """获取所有插件"""
        return list(self._plugins)


_global_registry: Optional[PluginRegistry] = None
_global_loader: Optional[PluginLoader] = None


def get_plugin_registry() -> PluginRegistry:
    """获取全局插件注册表"""
    global _global_registry
    if _global_registry is None:
        _global_registry = PluginRegistry()
    return _global_registry


def get_plugin_loader() -> PluginLoader:
    """获取全局插件加载器"""
    global _global_loader
    if _global_loader is None:
        _global_loader = PluginLoader(get_plugin_registry())
    return _global_loader


def load_all_plugins(search_paths: Optional[List[Path]] = None) -> int:
    """加载所有可用插件"""
    loader = get_plugin_loader()
    return loader.load_all(search_paths)


def register_qc_plugin(plugin: QCPlugin) -> None:
    """注册质控插件到全局注册表"""
    registry = get_plugin_registry()
    metadata = plugin.get_metadata()
    registry.register_plugin(
        metadata.name,
        LoadedPlugin(metadata=metadata, instance=plugin, module_path="<builtin>"),
    )


# ==================== 分类学解析插件钩子 ====================

_taxonomy_parsers: List[Callable] = []


def register_taxonomy_parser(parser_func: Callable) -> None:
    """注册自定义分类学解析器插件。

    Args:
        parser_func: 接受 (label: str, mode: str) -> dict 的函数
    """
    _taxonomy_parsers.append(parser_func)
    # 同时注册到 QC 插件链（如果需要）
    logger.info(f"Registered taxonomy parser plugin: {parser_func.__name__}")


def get_taxonomy_parsers() -> List[Callable]:
    """获取所有已注册的分类学解析插件"""
    return list(_taxonomy_parsers)


# ==================== 输出格式生成插件钩子 ====================

_output_formatters: List[Callable] = []


def register_output_formatter(formatter_func: Callable) -> None:
    """注册自定义输出格式生成器插件。

    Args:
        formatter_func: 接受 (results: dict, output_path: Path) -> None 的函数
    """
    _output_formatters.append(formatter_func)
    logger.info(f"Registered output formatter plugin: {formatter_func.__name__}")


def get_output_formatters() -> List[Callable]:
    """获取所有已注册的输出格式生成器插件"""
    return list(_output_formatters)


# ==================== 树验证插件钩子 ====================

_tree_validators: List[Callable] = []


def register_tree_validator(validator_func: Callable) -> None:
    """注册自定义树验证插件。

    Args:
        validator_func: 接受 (tree_str: str) -> List[str] 的函数，返回错误列表
    """
    _tree_validators.append(validator_func)
    logger.info(f"Registered tree validator plugin: {validator_func.__name__}")


def get_tree_validators() -> List[Callable]:
    """获取所有已注册的树验证插件"""
    return list(_tree_validators)


def create_example_qc_plugin() -> str:
    """生成示例 QC 插件代码"""
    return '''"""
Long Branch Attraction QC Plugin

移除可能导致长枝吸引假说的序列
"""

from phylodater.infrastructure.plugins import (
    QCPlugin, PluginMetadata, PluginType
)

PLUGIN_METADATA = {
    "name": "lba_filter",
    "version": "0.1.0",
    "description": "Remove sequences that may cause Long Branch Attraction artifacts",
    "plugin_type": "qc",
    "author": "Example Author",
    "dependencies": [],
    "config_schema": {
        "type": "object",
        "properties": {
            "threshold": {
                "type": "number",
                "description": "Branch length threshold for LBA detection",
                "default": 1.0
            }
        }
    },
    "config": {}
}

class PluginClass(QCPlugin):
    def get_metadata(self):
        return PluginMetadata(**PLUGIN_METADATA)

    def initialize(self, config):
        self.threshold = config.get("threshold", 1.0)

    def process_tree(self, tree, context):
        """移除长枝"""
        # 实现 LBA 检测和移除逻辑
        return tree

    def process_alignment(self, alignment, context):
        """移除可能导致 LBA 的序列"""
        # 实现序列移除逻辑
        return alignment
'''
