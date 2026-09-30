"""
TreeValidator 服务单元测试

测试树结构验证服务的各项功能
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock

import pytest

from phylodater.models import PhylogeneticTree
from phylodater.services.tree_validator import TreeValidator, ValidationReport


class TestTreeValidator:
    """TreeValidator 测试类"""

    @pytest.fixture
    def validator(self):
        """创建验证器实例"""
        return TreeValidator()

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        return tree

    def test_validate_rooted_true(self, validator, mock_tree):
        """测试定根验证 - 已定时"""
        mock_tree.is_rooted.return_value = True

        result = validator.validate_rooted(mock_tree)

        assert result is True

    def test_validate_rooted_false(self, validator, mock_tree):
        """测试定根验证 - 未定时"""
        mock_tree.is_rooted.return_value = False

        result = validator.validate_rooted(mock_tree)

        assert result is False

    @pytest.mark.parametrize(
        "newick,expected",
        [
            ("((A:0.1,B:0.1):0.2,(C:0.3,D:0.3):0.1);", True),
            ("((A:0.1,B:0.1):0.2,(C:0.3,D:0.3):0.1,(E:0.2,F:0.2):0.1);", False),
        ],
    )
    def test_validate_rooted_with_real_tree(self, validator, newick, expected):
        """真实 ``PhylogeneticTree.is_rooted`` 是**属性**而不是方法。

        旧实现写 ``tree.is_rooted()`` → 对真树必然 ``TypeError: 'bool' object is
        not callable``（只有 Mock 才让它看起来能用）。
        """
        tree = PhylogeneticTree(newick)

        assert validator.validate_rooted(tree) is expected

    def test_validate_binary_true(self, validator, mock_tree):
        """测试二叉性验证 - 二叉树"""
        mock_tree.is_binary.return_value = True

        result = validator.validate_binary(mock_tree)

        assert result is True

    def test_validate_binary_false(self, validator, mock_tree):
        """测试二叉性验证 - 非二叉树"""
        mock_tree.is_binary.return_value = False

        result = validator.validate_binary(mock_tree)

        assert result is False

    def test_validate_branch_lengths_valid(self, validator, mock_tree):
        """测试分支长度验证 - 有效"""
        mock_tree.validate_branch_lengths.return_value = (True, [], [])

        is_valid, errors, warnings = validator.validate_branch_lengths(mock_tree)

        assert is_valid is True
        assert errors == []
        assert warnings == []

    def test_validate_branch_lengths_invalid(self, validator, mock_tree):
        """测试分支长度验证 - 无效"""
        mock_tree.validate_branch_lengths.return_value = (
            False,
            ["Negative branch length found"],
            [],
        )

        is_valid, errors, warnings = validator.validate_branch_lengths(mock_tree)

        assert is_valid is False
        assert len(errors) == 1

    def test_normalize_name(self, validator):
        """测试名称规范化"""
        assert validator._normalize_name("Human_1") == "human_1"
        assert validator._normalize_name("Seq;Name") == "seq_name"
        assert validator._normalize_name("  Test  ") == "test"

    def test_full_validation_pass(self, validator, mock_tree):
        """测试完整验证 - 通过"""
        mock_tree.get_rooting_info.return_value = {
            "is_rooted": True,
            "likely_artifact": False,
        }
        mock_tree.is_binary.return_value = True
        mock_tree.validate_branch_lengths.return_value = (True, [], [])

        report = validator.full_validation(mock_tree)

        assert report.is_valid is True
        assert len(report.errors) == 0

    def test_full_validation_fail_unrooted(self, validator, mock_tree):
        """测试完整验证 - 未定时失败"""
        mock_tree.get_rooting_info.return_value = {
            "is_rooted": False,
            "likely_artifact": False,
        }
        mock_tree.is_binary.return_value = True
        mock_tree.validate_branch_lengths.return_value = (True, [], [])

        report = validator.full_validation(mock_tree, require_rooted=True)

        assert report.is_valid is False
        assert any("not rooted" in e.lower() for e in report.errors)

    def test_full_validation_with_alignment(self, validator, mock_tree):
        """测试带比对文件的完整验证"""
        mock_tree.get_rooting_info.return_value = {
            "is_rooted": True,
            "likely_artifact": False,
        }
        mock_tree.is_binary.return_value = True
        mock_tree.validate_branch_lengths.return_value = (True, [], [])

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">human\nATCG\n>chimp\nATCG\n>mouse\nATCG\n>rat\nATCG\n")
            alignment_path = Path(f.name)

        try:
            report = validator.full_validation(mock_tree, alignment_path)

            assert isinstance(report, ValidationReport)
            assert report.is_valid is True
        finally:
            alignment_path.unlink()


class TestSequenceConsistency:
    """序列一致性验证测试类"""

    @pytest.fixture
    def validator(self):
        return TreeValidator()

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        return tree

    def test_consistency_pass_fasta(self, validator, mock_tree):
        """测试 FASTA 格式一致性验证 - 通过"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">human\nATCG\n>chimp\nATCG\n>mouse\nATCG\n>rat\nATCG\n")
            alignment_path = Path(f.name)

        try:
            report = validator.validate_sequence_consistency(alignment_path, mock_tree)
            assert report.is_valid is True
            assert len(report.errors) == 0
        finally:
            alignment_path.unlink()

    def test_consistency_fail_count_mismatch(self, validator, mock_tree):
        """测试物种数量不一致"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">human\nATCG\n>chimp\nATCG\n>mouse\nATCG\n")
            alignment_path = Path(f.name)

        try:
            report = validator.validate_sequence_consistency(alignment_path, mock_tree)
            assert report.is_valid is False
            assert any("Species count mismatch" in e for e in report.errors)
            assert any("4 tips" in e for e in report.errors)
            assert any("3 sequences" in e for e in report.errors)
        finally:
            alignment_path.unlink()

    def test_consistency_fail_tree_extra(self, validator, mock_tree):
        """测试树中有比对中没有的物种"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">human\nATCG\n>chimp\nATCG\n>mouse\nATCG\n>dog\nATCG\n")
            alignment_path = Path(f.name)

        try:
            report = validator.validate_sequence_consistency(alignment_path, mock_tree)
            assert report.is_valid is False
            assert any(
                "tip(s) in tree not found in alignment" in e for e in report.errors
            )
            assert any("- rat" in e for e in report.errors)
        finally:
            alignment_path.unlink()

    def test_consistency_fail_aln_extra(self, validator, mock_tree):
        """测试比对中有树中没有的序列"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            # 比对包含 human, chimp, dog, cat；树包含 human, chimp, mouse, rat
            f.write(">human\nATCG\n>chimp\nATCG\n>dog\nATCG\n>cat\nATCG\n")
            alignment_path = Path(f.name)

        try:
            report = validator.validate_sequence_consistency(alignment_path, mock_tree)
            assert report.is_valid is False
            # 比对中有 dog, cat 不在树中
            assert any(
                "sequence(s) in alignment not found in tree" in e for e in report.errors
            )
            assert any("- cat" in e for e in report.errors)
            assert any("- dog" in e for e in report.errors)
            # 树中有 mouse, rat 不在比对中
            assert any(
                "tip(s) in tree not found in alignment" in e for e in report.errors
            )
            assert any("- mouse" in e for e in report.errors)
            assert any("- rat" in e for e in report.errors)
        finally:
            alignment_path.unlink()

    def test_consistency_fail_both_directions(self, validator, mock_tree):
        """测试双向不匹配"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">human\nATCG\n>chimp\nATCG\n>dog\nATCG\n>cat\nATCG\n")
            alignment_path = Path(f.name)

        try:
            report = validator.validate_sequence_consistency(alignment_path, mock_tree)
            assert report.is_valid is False
            # 应该有两个方向的错误
            tree_extra_errors = [
                e for e in report.errors if "tip(s) in tree not found" in e
            ]
            aln_extra_errors = [
                e for e in report.errors if "sequence(s) in alignment not found" in e
            ]
            assert len(tree_extra_errors) > 0
            assert len(aln_extra_errors) > 0
        finally:
            alignment_path.unlink()

    def test_consistency_with_phylip(self, validator, mock_tree):
        """测试 PHYLIP 格式一致性验证"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".phylip", delete=False) as f:
            # PHYLIP 格式：名称和序列之间用空格分隔
            f.write("4 4\nhuman ATCG\nchimp ATCG\nmouse ATCG\nrat   ATCG\n")
            alignment_path = Path(f.name)

        try:
            report = validator.validate_sequence_consistency(alignment_path, mock_tree)
            assert report.is_valid is True
        finally:
            alignment_path.unlink()

    def test_detect_alignment_format_fasta(self, validator):
        """测试 FASTA 格式检测"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\nATCG\n")
            path = Path(f.name)

        try:
            fmt = validator._detect_alignment_format(path)
            assert fmt == "fasta"
        finally:
            path.unlink()

    def test_detect_alignment_format_phylip(self, validator):
        """测试 PHYLIP 格式检测"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".phy", delete=False) as f:
            f.write("2 4\nseq1  ATCG\nseq2  ATCG\n")
            path = Path(f.name)

        try:
            fmt = validator._detect_alignment_format(path)
            assert fmt == "phylip"
        finally:
            path.unlink()

    def test_normalize_name(self, validator):
        """测试名称规范化在一致性检查中的应用"""
        # 带特殊字符的名称应该被规范化
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">Human_1\nATCG\n>Chimp;2\nATCG\n>Mouse(3)\nATCG\n>Rat:4\nATCG\n")
            alignment_path = Path(f.name)

        tree = Mock(spec=PhylogeneticTree)
        tree.tip_names = ["Human_1", "Chimp;2", "Mouse(3)", "Rat:4"]

        try:
            report = validator.validate_sequence_consistency(alignment_path, tree)
            assert report.is_valid is True
        finally:
            alignment_path.unlink()


class TestValidationReport:
    """ValidationReport 测试类"""

    def test_valid_report(self):
        """测试有效报告"""
        report = ValidationReport(is_valid=True, errors=[], warnings=[])

        assert report.is_valid is True
        assert len(report.errors) == 0
        assert len(report.warnings) == 0


class TestTimeTreeAssumption:
    """审阅项 B-18：输入树是否为时间树（超度量），按方法分档报告。"""

    ULTRAMETRIC = "((A:0.5,B:0.5):0.3,(C:0.4,D:0.4):0.4);"
    NON_ULTRAMETRIC = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"

    @pytest.fixture
    def validator(self):
        return TreeValidator()

    def _tree(self, newick):
        return PhylogeneticTree(newick)

    def test_ultrametric_input_produces_no_fuss(self, validator):
        report = validator.validate_ultrametric(
            self._tree(self.ULTRAMETRIC), method="lsd2"
        )

        assert report.is_valid is True
        assert report.errors == [] and report.warnings == []

    @pytest.mark.parametrize("method", ["lsd2", "iqtree2", "mcmctree", None])
    def test_non_ultrametric_warns_for_time_tree_assuming_methods(
        self, validator, method
    ):
        """LSD2 / 尖端定年 / MCMCTree 全局钟：warn（但不 hard-fail）。"""
        report = validator.validate_ultrametric(
            self._tree(self.NON_ULTRAMETRIC), method=method
        )

        assert report.is_valid is True
        assert len(report.warnings) == 1
        assert "超度量" in report.warnings[0] or "时间树" in report.warnings[0]

    @pytest.mark.parametrize("method", ["treepl", "pathd8", "r8s", "pyr8s"])
    def test_non_ultrametric_is_info_only_for_tolerant_methods(self, validator, method):
        """treePL / PATHd8 / r8s 设计上就吃非超度量树 → 不进 warning 通道。"""
        report = validator.validate_ultrametric(
            self._tree(self.NON_ULTRAMETRIC), method=method
        )

        assert report.is_valid is True
        assert report.warnings == []
        assert report.errors == []

    def test_require_time_tree_escalates_to_error(self, validator):
        report = validator.validate_ultrametric(
            self._tree(self.NON_ULTRAMETRIC), method="lsd2", require_time_tree=True
        )

        assert report.is_valid is False
        assert len(report.errors) == 1

    def test_full_validation_adds_time_tree_warning_without_failing(self, validator):
        """流水线 ``_validate_inputs`` 会在 ``is_valid=False`` 时抛错——
        所以 B-18 这道新检查默认只能进 warning 通道，否则所有非超度量输入
        都会被硬拦（对 treePL/PATHd8 用户是回归）。"""
        report = validator.full_validation(self._tree(self.NON_ULTRAMETRIC))

        assert report.is_valid is True
        assert any("时间树" in w or "超度量" in w for w in report.warnings)

    def test_full_validation_can_hard_require_time_tree(self, validator):
        report = validator.full_validation(
            self._tree(self.NON_ULTRAMETRIC),
            method="lsd2",
            require_time_tree=True,
        )

        assert report.is_valid is False
        assert any("超度量" in e or "时间树" in e for e in report.errors)

    def test_full_validation_passes_quietly_on_time_tree(self, validator):
        report = validator.full_validation(self._tree(self.ULTRAMETRIC), method="lsd2")

        assert report.is_valid is True
        assert not any("时间树" in w for w in report.warnings)

    def test_cladogram_is_reported_as_undecidable(self, validator):
        """没有分支长度的拓扑树：既不能说"是时间树"，也不该判死。"""
        report = validator.validate_ultrametric(
            self._tree("((A,B),(C,D));"), method="lsd2"
        )

        assert report.is_valid is True
        assert len(report.warnings) == 1
        assert "无法判定" in report.warnings[0]


class TestValidationReportShapes:
    """ValidationReport 数据类本身的形状约定。"""

    def test_invalid_report_with_errors(self):
        """测试带错误的无效报告"""
        report = ValidationReport(
            is_valid=False, errors=["Error 1", "Error 2"], warnings=["Warning 1"]
        )

        assert report.is_valid is False
        assert len(report.errors) == 2
        assert len(report.warnings) == 1

    def test_report_with_only_warnings(self):
        """测试仅带警告的报告"""
        report = ValidationReport(
            is_valid=True, errors=[], warnings=["This is a warning"]
        )

        assert report.is_valid is True
        assert len(report.warnings) == 1
