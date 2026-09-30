"""
模块二：输入文件验证功能测试

验证树文件、序列文件、交叉验证、编码兼容等场景。
"""

import pytest

from phylodater import PhylogeneticTree
from phylodater.services.deep_validator import DeepValidator, SequenceAlphabet
from phylodater.services.input_validator import InputValidator


class TestTreeFileValidation:
    """树文件验证测试"""

    def test_val_001_valid_newick(self, normal_tree):
        """VAL-001: 有效 Newick"""
        tree = PhylogeneticTree.from_file(str(normal_tree))
        assert tree is not None
        assert tree.num_tips > 0

    def test_val_002_valid_nexus(self, functional_input_dir):
        """VAL-002: 有效 Nexus"""
        nexus_file = functional_input_dir / "boundary" / "multiple_trees.nex"
        tree = PhylogeneticTree.from_file(str(nexus_file), multi_tree_mode="first")
        assert tree is not None

    def test_val_003_file_not_found(self):
        """VAL-003: 文件不存在"""
        with pytest.raises(FileNotFoundError):
            PhylogeneticTree.from_file("/nonexistent/tree.nwk")

    def test_val_004_empty_file(self, functional_input_dir):
        """VAL-004: 空文件"""
        empty_tree = functional_input_dir / "boundary" / "empty_tree.nwk"
        validator = DeepValidator()
        result = validator.validate_tree_deep(empty_tree)
        assert not result.is_valid

    def test_val_005_unknown_format(self, tmp_path):
        """VAL-005: 未知格式"""
        bad_file = tmp_path / "random.txt"
        bad_file.write_text("this is not a tree", encoding="utf-8")
        validator = InputValidator()
        result = validator.validate_tree_file(bad_file)
        assert not result.is_valid or "无法检测" in str(result.errors)

    def test_val_010_unbalanced_parentheses(self, functional_input_dir):
        """VAL-010: 括号不平衡"""
        unbalanced = functional_input_dir / "invalid" / "unbalanced_tree.nwk"
        validator = DeepValidator()
        result = validator.validate_tree_deep(unbalanced)
        assert not result.is_valid
        assert any("括号" in e for e in result.errors)

    def test_val_011_extra_parentheses(self, tmp_path):
        """VAL-011: 多余右括号"""
        bad_tree = tmp_path / "extra.nwk"
        bad_tree.write_text("((A:0.1,B:0.2):0.3,C:0.4));", encoding="utf-8")
        validator = DeepValidator()
        result = validator.validate_tree_deep(bad_tree)
        assert not result.is_valid

    def test_val_012_negative_branch_length(self, functional_input_dir):
        """VAL-012: 负分支长度应抛出异常"""
        from phylodater.core.exceptions import NegativeBranchLengthError

        negative = functional_input_dir / "invalid" / "negative_branch_tree.nwk"
        validator = DeepValidator()
        with pytest.raises(NegativeBranchLengthError):
            validator.validate_tree_deep(negative)

    def test_val_013_duplicate_leaf_names(self, functional_input_dir):
        """VAL-013: 重复末端节点名"""
        duplicate = functional_input_dir / "invalid" / "duplicate_leaf_names.nwk"
        validator = DeepValidator()
        result = validator.validate_tree_deep(duplicate)
        assert any(
            "重复" in e or "duplicate" in e.lower() for e in result.errors
        ) or any("重复" in w or "duplicate" in w.lower() for w in result.warnings)

    def test_val_014_duplicate_internal_labels(self, functional_input_dir):
        """VAL-014: 重复内部节点名"""
        duplicate_internal = (
            functional_input_dir / "invalid" / "duplicate_internal_labels.nwk"
        )
        validator = DeepValidator()
        result = validator.validate_tree_deep(duplicate_internal)
        assert any(
            "重复" in w or "duplicate" in w.lower() for w in result.warnings
        ) or any("重复" in e or "duplicate" in e.lower() for e in result.errors)

    def test_val_015_empty_node_name(self, tmp_path):
        """VAL-015: 空节点名"""
        bad_tree = tmp_path / "empty_node.nwk"
        bad_tree.write_text("((:0.1,B:0.2):0.3,C:0.4);", encoding="utf-8")
        validator = DeepValidator()
        result = validator.validate_tree_deep(bad_tree)
        assert not result.is_valid or any(
            "空" in e or "empty" in e.lower() for e in result.errors
        )

    @pytest.mark.parametrize(
        "filename,desc",
        [
            ("control_char_tree.nwk", "控制字符"),
            ("unicode_bidi_tree.nwk", "Unicode 双向文本"),
            ("zero_width_tree.nwk", "零宽字符"),
        ],
    )
    def test_val_016_017_malicious_chars(self, functional_input_dir, filename, desc):
        """VAL-016/017: 恶意字符检测"""
        bad_tree = functional_input_dir / "invalid" / filename
        validator = DeepValidator()
        result = validator.validate_tree_deep(bad_tree)
        assert not result.is_valid or any(
            "恶意" in e or "invalid" in e.lower() for e in result.errors
        )


class TestSequenceFileValidation:
    """序列文件验证测试"""

    def test_val_020_valid_fasta(self, normal_alignment):
        """VAL-020: 有效 FASTA"""
        from Bio import SeqIO

        records = list(SeqIO.parse(str(normal_alignment), "fasta"))
        assert len(records) == 5

    def test_val_022_empty_alignment(self, functional_input_dir):
        """VAL-022: 空序列文件"""
        empty_aln = functional_input_dir / "invalid" / "empty_alignment.fasta"
        validator = DeepValidator()
        result = validator.validate_alignment_deep(empty_aln)
        assert not result.is_valid or result.num_sequences == 0

    def test_val_023_duplicate_seq_id(self, functional_input_dir):
        """VAL-023: 重复序列 ID"""
        dup_aln = functional_input_dir / "invalid" / "duplicate_seq_id.fasta"
        validator = DeepValidator()
        result = validator.validate_alignment_deep(dup_aln)
        assert len(result.duplicate_ids) > 0

    def test_val_024_invalid_dna_chars(self, functional_input_dir):
        """VAL-024: DNA 无效字符"""
        invalid_aln = functional_input_dir / "invalid" / "invalid_dna_chars.fasta"
        validator = DeepValidator()
        result = validator.validate_alignment_deep(
            invalid_aln, expected_alphabet=SequenceAlphabet.DNA
        )
        assert not result.is_valid or len(result.invalid_chars) > 0

    def test_val_026_phylip_format(self, functional_input_dir):
        """VAL-026: PHYLIP 格式检测"""
        phylip = functional_input_dir / "invalid" / "phylip_format.phy"
        validator = InputValidator()
        result = validator.validate_alignment_file(phylip)
        # 当前代码可能返回 phylip 有效或提示不支持
        assert result.format_detected == "phylip" or any(
            "不支持" in e for e in result.errors
        )

    def test_val_027_stockholm_format(self, functional_input_dir):
        """VAL-027: Stockholm 格式检测"""
        sto = functional_input_dir / "invalid" / "stockholm_format.sto"
        validator = InputValidator()
        result = validator.validate_alignment_file(sto)
        assert result.format_detected == "stockholm" or any(
            "不支持" in e for e in result.errors
        )


class TestCrossValidation:
    """交叉验证测试"""

    def test_val_030_label_match(self, normal_tree, normal_alignment):
        """VAL-030: 标签完全匹配"""
        from phylodater.services.tree_validator import TreeValidator

        tree = PhylogeneticTree.from_file(str(normal_tree))
        validator = TreeValidator()
        report = validator.validate_sequence_consistency(normal_alignment, tree)
        assert report.is_valid

    def test_val_031_tree_extra_label(self, tmp_path):
        """VAL-031: 树中多出标签"""
        from phylodater.services.tree_validator import TreeValidator

        tree = PhylogeneticTree.from_newick("((A:0.1,B:0.2):0.3,C:0.4);")
        aln = tmp_path / "aln.fasta"
        aln.write_text(">A\nATCG\n>B\nATCG\n", encoding="utf-8")
        validator = TreeValidator()
        report = validator.validate_sequence_consistency(aln, tree)
        assert not report.is_valid
        assert any(
            "tree" in e.lower() and "not found" in e.lower() for e in report.errors
        )

    def test_val_032_alignment_extra_label(self, tmp_path):
        """VAL-032: 序列中多出标签"""
        from phylodater.services.tree_validator import TreeValidator

        tree = PhylogeneticTree.from_newick("((A:0.1,B:0.2):0.3,C:0.4);")
        aln = tmp_path / "aln.fasta"
        aln.write_text(">A\nATCG\n>B\nATCG\n>D\nATCG\n", encoding="utf-8")
        validator = TreeValidator()
        report = validator.validate_sequence_consistency(aln, tree)
        assert not report.is_valid
        assert any(
            "alignment" in e.lower() or "sequence" in e.lower() for e in report.errors
        )


class TestEncodingCompatibility:
    """编码兼容性测试"""

    def test_val_040_gbk_encoding(self, functional_input_dir):
        """VAL-040: GBK 编码"""
        gbk_tree = functional_input_dir / "boundary" / "gbk_tree.nwk"
        tree = PhylogeneticTree.from_file(str(gbk_tree))
        assert tree is not None

    def test_val_041_bom_header(self, functional_input_dir):
        """VAL-041: BOM 头"""
        bom_tree = functional_input_dir / "boundary" / "bom_tree.nwk"
        tree = PhylogeneticTree.from_file(str(bom_tree))
        assert tree is not None

    def test_val_042_chinese_path(self, tmp_path):
        """VAL-042: 中文路径"""

        chinese_dir = tmp_path / "中文路径"
        chinese_dir.mkdir()
        tree_file = chinese_dir / "tree.nwk"
        tree_file.write_text("((A:0.1,B:0.2):0.3,C:0.4);", encoding="utf-8")
        tree = PhylogeneticTree.from_file(str(tree_file))
        assert tree is not None

    def test_val_043_space_path(self, tmp_path):
        """VAL-043: 空格路径"""
        space_dir = tmp_path / "space path"
        space_dir.mkdir()
        tree_file = space_dir / "tree.nwk"
        tree_file.write_text("((A:0.1,B:0.2):0.3,C:0.4);", encoding="utf-8")
        tree = PhylogeneticTree.from_file(str(tree_file))
        assert tree is not None
