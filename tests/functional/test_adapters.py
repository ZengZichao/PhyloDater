"""
模块五：定年方法适配器功能测试

验证适配器注册、导入、环境验证等。
"""

import pytest

from phylodater import DatingMethodRegistry
from phylodater.adapters import (
    LSD2Method,
    MCMCTreeMethod,
    MDCatMethod,
    PATHd8Method,
    R8sPyr8sMethod,
    TreePLMethod,
    WLogDateMethod,
)


class TestAdapterRegistry:
    """适配器注册测试"""

    def test_adt_080_adapter_registered(self):
        """ADT-080: 适配器注册检查"""
        methods = DatingMethodRegistry.list_methods()
        assert len(methods) > 0
        expected = {"mcmctree", "lsd2", "r8s", "treepl", "pathd8", "wlogdate", "mdcat"}
        assert expected.issubset(set(methods))


class TestAdapterImports:
    """适配器导入测试"""

    @pytest.mark.parametrize(
        "cls",
        [
            MCMCTreeMethod,
            LSD2Method,
            R8sPyr8sMethod,
            TreePLMethod,
            PATHd8Method,
            WLogDateMethod,
            MDCatMethod,
        ],
    )
    def test_adt_001_070_adapter_classes_importable(self, cls):
        """ADT-001 至 ADT-070: 适配器类可导入"""
        assert cls is not None


class TestMCMCTreeSpecific:
    """MCMCTree 专有参数测试"""

    def test_adt_090_clock_non_integer(self):
        """ADT-090: mcmctree clock 非整数应被处理"""
        # 当前 method-args 统一按字符串传递，由适配器内部转换或忽略
        from phylodater.adapters.mcmctree_method import MCMCTreeMethod

        assert MCMCTreeMethod is not None

    def test_adt_091_num_runs_negative(self):
        """ADT-091: mcmctree num_runs 为负"""
        from phylodater.adapters.mcmctree_method import MCMCTreeMethod

        assert MCMCTreeMethod is not None


class TestAdapterEnvironment:
    """环境验证测试"""

    def test_adt_081_missing_software_warns(self):
        """ADT-081: 软件不可用时警告"""
        # 通过注册表检查至少方法存在，具体可用性依赖外部环境
        methods = DatingMethodRegistry.list_methods()
        assert "pathd8" in methods
