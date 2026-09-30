"""
模块八：自检模式功能测试

验证依赖检查、功能测试、输出格式。
"""

from tests.functional.conftest import run_phylodater_cli


class TestSelfCheck:
    """自检模式测试"""

    def test_chk_001_010_dependencies(self):
        """CHK-001 至 CHK-010: 依赖检查"""
        result = run_phylodater_cli(["check"])
        assert result.returncode in [0, 1]
        assert "[PASS]" in result.stdout or "[FAIL]" in result.stdout

    def test_chk_020_023_functional_tests(self):
        """CHK-020 至 CHK-023: 功能测试"""
        result = run_phylodater_cli(["check"])
        assert result.returncode in [0, 1]
        assert (
            "树解析" in result.stdout
            or "分类学" in result.stdout
            or "单系群" in result.stdout
        )

    def test_chk_040_table_format(self):
        """CHK-040/042: 表格格式输出"""
        result = run_phylodater_cli(["check"])
        assert result.returncode in [0, 1]
        assert "[PASS]" in result.stdout or "[FAIL]" in result.stdout
        assert "总计:" in result.stdout
