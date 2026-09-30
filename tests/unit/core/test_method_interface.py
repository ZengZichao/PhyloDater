"""
DatingMethod 接口和 Registry 的单元测试
"""

import dendropy
import pytest

from phylodater.core.exceptions import UnknownMethodError
from phylodater.core.method_interface import DatingMethod, DatingMethodRegistry
from phylodater.infrastructure.configuration import (
    CommonConfig,
    SoftwarePaths,
    ToolConfig,
)
from phylodater.models.calibration import CalibrationPoint
from phylodater.models.constraints import FixedAgeConstraint


class MockDatingMethod(DatingMethod):
    """模拟的定年方法实现"""

    @property
    def method_name(self) -> str:
        return "mock_method"

    def validate_environment(self) -> bool:
        return True

    def prepare_inputs(self, tree, calibrations, alignment_path=None):
        return {"prepared": True}

    def execute(self) -> bool:
        return True

    def parse_results(self):
        from phylodater.models import DatingResult

        return DatingResult(method_name="mock_method", node_ages={})


class MockDatingMethodWithCommonConfig(DatingMethod):
    """接受 common_config 的模拟定年方法实现"""

    def __init__(self, config, output_dir, software_paths=None, common_config=None):
        super().__init__(config, output_dir, software_paths)
        self.common_config = common_config or CommonConfig()

    @property
    def method_name(self) -> str:
        return "mock_method_common"

    def validate_environment(self) -> bool:
        return True

    def prepare_inputs(self, tree, calibrations, alignment_path=None):
        return {"prepared": True}

    def execute(self) -> bool:
        return True

    def parse_results(self):
        from phylodater.models import DatingResult

        return DatingResult(method_name="mock_method_common", node_ages={})


class TestDatingMethodRegistry:
    """测试 DatingMethodRegistry"""

    def setup_method(self):
        """每个测试前清理注册表"""
        DatingMethodRegistry._registry.clear()

    def teardown_method(self):
        """每个测试后清理注册表"""
        DatingMethodRegistry._registry.clear()

    def test_register_method(self):
        """测试注册方法"""
        DatingMethodRegistry.register("mock", MockDatingMethod)
        assert "mock" in DatingMethodRegistry.list_methods()

    def test_register_case_insensitive(self):
        """测试注册不区分大小写"""
        DatingMethodRegistry.register("MockMethod", MockDatingMethod)
        assert "mockmethod" in DatingMethodRegistry.list_methods()

    def test_create_method(self, tmp_path):
        """测试创建方法实例"""
        DatingMethodRegistry.register("mock", MockDatingMethod)

        config = ToolConfig()
        # tmp_path，不是硬编码 /tmp：适配器会在 output_dir 下创建工作目录，
        # 写死的共享路径既依赖 /tmp 可写，又会与并行作业相互覆盖。
        output_dir = tmp_path

        instance = DatingMethodRegistry.create("mock", config, output_dir)

        assert isinstance(instance, MockDatingMethod)
        assert instance.method_name == "mock_method"

    def test_create_method_with_common_config(self, tmp_path):
        """测试创建方法实例并传递 common_config"""
        DatingMethodRegistry.register("mock_common", MockDatingMethodWithCommonConfig)

        config = ToolConfig()
        config.common.nthreads = 8
        config.common.seed = 42
        output_dir = tmp_path

        instance = DatingMethodRegistry.create("mock_common", config, output_dir)

        assert isinstance(instance, MockDatingMethodWithCommonConfig)
        assert instance.common_config.nthreads == 8
        assert instance.common_config.seed == 42

    def test_create_unknown_method_raises(self, tmp_path):
        """测试创建未知方法应抛出异常"""
        config = ToolConfig()
        output_dir = tmp_path

        with pytest.raises(UnknownMethodError, match="Unknown dating method"):
            DatingMethodRegistry.create("unknown", config, output_dir)

    def test_is_registered(self):
        """测试检查方法是否已注册"""
        DatingMethodRegistry.register("mock", MockDatingMethod)

        assert DatingMethodRegistry.is_registered("mock") is True
        assert DatingMethodRegistry.is_registered("MOCK") is True  # 不区分大小写
        assert DatingMethodRegistry.is_registered("unknown") is False

    def test_list_methods(self):
        """测试列出所有注册的方法"""
        DatingMethodRegistry.register("method1", MockDatingMethod)
        DatingMethodRegistry.register("method2", MockDatingMethod)

        methods = DatingMethodRegistry.list_methods()
        assert "method1" in methods
        assert "method2" in methods
        assert len(methods) == 2


class TestDatingMethod:
    """测试 DatingMethod 基类"""

    def test_init_creates_work_dir(self, tmp_path):
        """测试初始化创建工作目录"""
        config = ToolConfig()
        output_dir = tmp_path / "output"

        method = MockDatingMethod(config, output_dir)

        assert method.work_dir.exists()
        assert method.work_dir.name.startswith("mock_method")

    def test_get_software_path(self, tmp_path):
        """测试获取软件路径"""
        config = ToolConfig()
        output_dir = tmp_path / "output"
        software_paths = SoftwarePaths(paml_path="/opt/paml")

        method = MockDatingMethod(config, output_dir, software_paths)

        assert method.get_software_path("paml_path") == "/opt/paml"
        assert method.get_software_path("nonexistent") is None

    def test_get_software_path_default(self, tmp_path):
        """测试默认软件路径为空"""
        config = ToolConfig()
        output_dir = tmp_path / "output"

        method = MockDatingMethod(config, output_dir)

        assert method.get_software_path("any_path") is None

    def test_extract_chronogram_ages_root_calibration(self, tmp_path):
        """根校准不应崩溃，且应返回正确的根年龄（≈3.0），不为 0。

        回归保护：dendropy 4.6.4 的 Node 没有 ``max_distance_from_tip``，
        旧写法会抛 AttributeError。改用 ``Tree.max_distance_from_root()`` 取根年龄。
        """
        config = ToolConfig()
        output_dir = tmp_path / "output"
        method = MockDatingMethod(config, output_dir)

        tree = dendropy.Tree.get(
            data="((A:1.0,B:1.0):2.0,(C:1.0,D:1.0):2.0);", schema="newick"
        )
        root_cal = CalibrationPoint(
            name="Root",
            is_root_node=True,
            age_constraint=FixedAgeConstraint(fixed_age=3.0),
        )

        # 不应抛 AttributeError
        node_ages = method._extract_chronogram_ages(tree, [root_cal], unit_factor=1.0)

        assert "Root" in node_ages
        assert node_ages["Root"].mean_age == pytest.approx(3.0, abs=1e-6)
        assert node_ages["Root"].mean_age != 0
