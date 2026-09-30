"""
DeepValidator 深度验证器单元测试

测试树文件和序列文件的深层验证功能
"""

import tempfile
from pathlib import Path

import pytest

from phylodater.services.deep_validator import DeepValidator, SequenceAlphabet


class TestDeepValidatorTreeValidation:
    """树文件深度验证测试"""

    @pytest.fixture
    def validator(self):
        return DeepValidator()

    @pytest.fixture
    def temp_dir(self):
        import shutil

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_valid_newick_tree(self, validator, temp_dir):
        """测试有效的 Newick 树"""
        tree_file = temp_dir / "valid.nwk"
        tree_file.write_text("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is True
        assert result.format_detected == "newick"
        assert result.tree_count == 1
        assert len(result.errors) == 0

    def test_bracket_imbalance_missing_closing(self, validator, temp_dir):
        """测试缺少右括号"""
        tree_file = temp_dir / "imbalance.nwk"
        tree_file.write_text("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6;")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is False
        assert any("缺少" in e and "右括号" in e for e in result.errors)

    def test_bracket_imbalance_extra_closing(self, validator, temp_dir):
        """测试多余的右括号"""
        tree_file = temp_dir / "extra.nwk"
        tree_file.write_text("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6));")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is False
        assert any("多余" in e and "右括号" in e for e in result.errors)

    def test_negative_branch_length(self, validator, temp_dir):
        """测试负分支长度 - 应抛出 NegativeBranchLengthError"""
        from phylodater.core.exceptions import NegativeBranchLengthError

        tree_file = temp_dir / "negative.nwk"
        tree_file.write_text("((A:0.1,B:-0.2):0.3,(C:0.4,D:0.5):0.6);")

        with pytest.raises(NegativeBranchLengthError):
            validator.validate_tree_deep(tree_file)

    def test_duplicate_node_names(self, validator, temp_dir):
        """测试重复节点名 - 末端重复应为ERROR"""
        tree_file = temp_dir / "duplicate.nwk"
        tree_file.write_text("((A:0.1,B:0.2):0.3,(A:0.4,D:0.5):0.6);")

        result = validator.validate_tree_deep(tree_file)

        # 重复的末端节点名应该是错误（会混淆分类学映射和序列匹配）
        assert any("重复的末端节点名" in e for e in result.errors)

    def test_multiple_trees(self, validator, temp_dir):
        """测试多棵树"""
        tree_file = temp_dir / "multi.nwk"
        tree_file.write_text(
            "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\n((E:0.1,F:0.2):0.3,(G:0.4,H:0.5):0.6);"
        )

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is True
        assert result.tree_count == 2
        assert len(result.tree_summaries) == 2
        assert any("2 棵树" in w for w in result.warnings)

    def test_nexus_format(self, validator, temp_dir):
        """测试 Nexus 格式"""
        tree_file = temp_dir / "test.nex"
        tree_file.write_text("""#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=4;
  TAXLABELS A B C D;
END;
BEGIN TREES;
  TREE tree1 = ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);
END;
""")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is True
        assert result.format_detected == "nexus"

    def test_empty_file(self, validator, temp_dir):
        """测试空文件"""
        tree_file = temp_dir / "empty.nwk"
        tree_file.write_text("")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is False
        assert any("为空" in e for e in result.errors)

    def test_nonexistent_file(self, validator, temp_dir):
        """测试不存在的文件"""
        tree_file = temp_dir / "nonexistent.nwk"

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is False
        assert any("不存在" in e for e in result.errors)

    def test_tree_summary_statistics(self, validator, temp_dir):
        """测试树的统计摘要"""
        tree_file = temp_dir / "summary.nwk"
        tree_file.write_text("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is True
        assert len(result.tree_summaries) == 1
        summary = result.tree_summaries[0]
        assert summary["terminals"] >= 2  # 至少2个末端节点
        assert summary["internal"] >= 1  # 至少1个内部节点


class TestDeepValidatorSequenceValidation:
    """序列文件深度验证测试"""

    @pytest.fixture
    def validator(self):
        return DeepValidator()

    @pytest.fixture
    def temp_dir(self):
        import shutil

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_valid_fasta(self, validator, temp_dir):
        """测试有效的 FASTA 文件"""
        fasta_file = temp_dir / "valid.fasta"
        fasta_file.write_text(">seq1\nATCGATCG\n>seq2\nATCGATCG\n")

        result = validator.validate_alignment_deep(fasta_file)

        assert result.is_valid is True
        assert result.format_detected == "fasta"
        assert result.num_sequences == 2
        assert result.sequence_length == 8
        assert result.is_aligned is True
        assert result.alphabet_detected == "dna"

    def test_duplicate_ids(self, validator, temp_dir):
        """测试重复序列ID"""
        fasta_file = temp_dir / "duplicate.fasta"
        fasta_file.write_text(">seq1\nATCG\n>seq1\nATCG\n")

        result = validator.validate_alignment_deep(fasta_file)

        assert result.is_valid is False
        assert len(result.duplicate_ids) == 1
        assert any("重复" in e for e in result.errors)

    def test_invalid_dna_chars(self, validator, temp_dir):
        """测试DNA序列中的无效字符"""
        fasta_file = temp_dir / "invalid.fasta"
        fasta_file.write_text(">seq1\nATCGXATC\n>seq2\nATCGATCG\n")

        result = validator.validate_alignment_deep(
            fasta_file, expected_alphabet=SequenceAlphabet.DNA
        )

        assert result.is_valid is False
        assert len(result.invalid_chars) > 0
        assert any("无效字符" in e for e in result.errors)

    def test_unequal_sequence_lengths(self, validator, temp_dir):
        """测试序列长度不一致"""
        fasta_file = temp_dir / "unequal.fasta"
        fasta_file.write_text(">seq1\nATCG\n>seq2\nATCGATCG\n")

        result = validator.validate_alignment_deep(fasta_file, require_aligned=True)

        assert result.is_valid is True  # 长度不一致是警告不是错误
        assert result.is_aligned is False
        assert any("长度不一致" in w for w in result.warnings)

    def test_protein_alphabet(self, validator, temp_dir):
        """测试蛋白质序列"""
        fasta_file = temp_dir / "protein.fasta"
        fasta_file.write_text(
            ">prot1\nMKWVTFISLLFLFSSAYS\n>prot2\nMKWVTFISLLFLFSSAYS\n"
        )

        result = validator.validate_alignment_deep(fasta_file)

        assert result.is_valid is True
        assert result.alphabet_detected == "protein"

    def test_rna_alphabet(self, validator, temp_dir):
        """测试RNA序列"""
        fasta_file = temp_dir / "rna.fasta"
        fasta_file.write_text(">rna1\nAUCGAUCG\n>rna2\nAUCGAUCG\n")

        result = validator.validate_alignment_deep(fasta_file)

        assert result.is_valid is True
        assert result.alphabet_detected == "rna"

    def test_fastq_format(self, validator, temp_dir):
        """测试FASTQ格式"""
        fastq_file = temp_dir / "test.fastq"
        fastq_file.write_text("@seq1\nATCG\n+\nIIII\n@seq2\nATCG\n+\nIIII\n")

        result = validator.validate_alignment_deep(fastq_file)

        assert result.is_valid is True
        assert result.format_detected == "fastq"
        assert result.num_sequences == 2

    def test_fastq_invalid_format(self, validator, temp_dir):
        """测试无效的FASTQ格式"""
        fastq_file = temp_dir / "invalid.fastq"
        fastq_file.write_text("@seq1\nATCG\n-\nIIII\n")  # 分隔符应该是+

        result = validator.validate_alignment_deep(fastq_file)

        assert result.is_valid is False

    def test_empty_fasta(self, validator, temp_dir):
        """测试空FASTA文件"""
        fasta_file = temp_dir / "empty.fasta"
        fasta_file.write_text("")

        result = validator.validate_alignment_deep(fasta_file)

        assert result.is_valid is False
        assert any("为空" in e for e in result.errors)

    def test_no_sequences(self, validator, temp_dir):
        """测试没有序列的FASTA文件"""
        fasta_file = temp_dir / "noseq.fasta"
        fasta_file.write_text(">header_only\n")

        result = validator.validate_alignment_deep(fasta_file)

        # 这应该被检测为空序列或无效
        assert result.num_sequences <= 1

    def test_nonexistent_file(self, validator, temp_dir):
        """测试不存在的文件"""
        fasta_file = temp_dir / "nonexistent.fasta"

        result = validator.validate_alignment_deep(fasta_file)

        assert result.is_valid is False
        assert any("不存在" in e for e in result.errors)

    def test_multiple_invalid_chars_with_line_numbers(self, validator, temp_dir):
        """测试多个无效字符的行号定位"""
        fasta_file = temp_dir / "multi_invalid.fasta"
        # 使用明确的无效字符 X, Y, Z（不是DNA字符）
        content = ">seq1\nATCGXATC\n>seq2\nATCGYZAT\n"
        fasta_file.write_text(content, encoding="utf-8")

        result = validator.validate_alignment_deep(
            fasta_file, expected_alphabet=SequenceAlphabet.DNA
        )

        # 如果有无效字符，验证行号信息
        if result.invalid_chars:
            for inv in result.invalid_chars:
                assert "line" in inv
                assert "char" in inv
                assert "id" in inv


class TestNodeNameValidity:
    """节点名合法性（wave-2a: P2 引号/数字开头标签）。

    Newick 规范允许：
    - 引号包裹的标签（可含空格/括号等），如 "sp A"
    - 以数字开头的标签，如 123sp
    仅未闭合引号等真正非法的情形才应被拒绝。
    """

    @pytest.fixture
    def validator(self):
        return DeepValidator()

    @pytest.fixture
    def temp_dir(self):
        import shutil

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_quoted_node_name_valid(self, validator, temp_dir):
        """带引号的节点名应通过校验（不卡死、不误判非法）。"""
        tree_file = temp_dir / "quoted.nwk"
        tree_file.write_text('((A:0.1,"sp B":0.2):0.3,(C:0.4,D:0.5):0.6);')
        result = validator.validate_tree_deep(tree_file)
        assert result.is_valid is True
        assert len(result.errors) == 0

    def test_digit_start_node_name_valid(self, validator, temp_dir):
        """以数字开头的节点名应通过校验。"""
        tree_file = temp_dir / "digit.nwk"
        tree_file.write_text("((123sp:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        result = validator.validate_tree_deep(tree_file)
        assert result.is_valid is True
        assert len(result.errors) == 0

    def test_quoted_name_with_leading_space_valid(self, validator, temp_dir):
        """引号标签前导空格不应被误判为空节点名。"""
        tree_file = temp_dir / "space.nwk"
        tree_file.write_text('(( "sp A":0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);')
        result = validator.validate_tree_deep(tree_file)
        assert result.is_valid is True
        assert len(result.errors) == 0

    def test_unclosed_quote_detected(self, validator, temp_dir):
        """未闭合的引号仍应被判定为非法（回归）。"""
        tree_file = temp_dir / "unclosed.nwk"
        tree_file.write_text('(( "sp A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);')
        result = validator.validate_tree_deep(tree_file)
        assert result.is_valid is False
        assert any("引号未正确关闭" in e for e in result.errors)

    def test_empty_node_name_detected(self, validator, temp_dir):
        """真正的空节点名（逗号间无标签）仍应被检测（回归）。"""
        tree_file = temp_dir / "empty.nwk"
        tree_file.write_text("((A:0.1,,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        result = validator.validate_tree_deep(tree_file)
        assert result.is_valid is False
        assert any("空节点名" in e for e in result.errors)


class TestBranchLengthsFromParsedTree:
    """审阅项 C-25 / C-26：分支长度检查必须读解析后的树，而不是扫原始文本。"""

    @pytest.fixture
    def validator(self):
        return DeepValidator()

    def test_negative_inside_quoted_leaf_name_is_not_a_branch_length(self, validator):
        """C-25 假阳性回归：引号叶名里的 ``:-3.5`` 不是分支长度。

        旧实现用正则 ``r":(-?\\d+\\.?\\d*...)"`` 直接扫文本，把完全合法的
        ``("A:-3.5":0.25,B:0.25);`` 报成 ``[CRITICAL] 负分支长度 -3.5``。
        """
        critical, warnings = validator._check_branch_lengths(
            '("A:-3.5":0.25,B:0.25)', "树"
        )

        assert critical == []
        assert warnings == []

    @pytest.mark.parametrize("newick", ["A:-.5,B:0.5", "A:-.5e1,B:0.5"])
    def test_leading_dot_negative_length_is_caught(self, validator, newick):
        """C-25 假阴性回归：``-.5`` 是 Newick 合法写法，旧正则完全不匹配。"""
        critical, _ = validator._check_branch_lengths(newick, "树")

        assert len(critical) == 1
        assert "CRITICAL" in critical[0]

    def test_negative_length_is_still_critical(self, validator):
        critical, warnings = validator._check_branch_lengths(
            "((A:0.1,B:-0.2):0.3,(C:0.4,D:0.5):0.6)", "树"
        )

        assert len(critical) == 1
        assert "-0.2" in critical[0]

    def test_zero_length_branch_populates_warning_channel(self, validator):
        """C-26：``warnings`` 通道此前恒为空，调用方以为拿到了一档软问题。"""
        critical, warnings = validator._check_branch_lengths("(A:0,B:0.5)", "树")

        assert critical == []
        assert any("零长度" in w for w in warnings)

    def test_absurd_length_populates_warning_channel(self, validator):
        """单位混淆（把 Ma 当替换/位点）最常见的表现是数量级异常，应留痕。"""
        critical, warnings = validator._check_branch_lengths("(A:2e9,B:0.5)", "树")

        assert critical == []
        assert any("异常大" in w for w in warnings)

    def test_unparseable_tree_degrades_loudly(self, validator):
        """解析不出来时才退回文本扫描，并且必须说明这是降级路径。"""
        critical, warnings = validator._check_branch_lengths("((A:0.1,B:0.2", "树")

        assert any("退回文本层扫描" in w for w in warnings)

    # ---- 同一件事走公开入口（validate_tree_deep）再钉一遍：报告实测的 FP/FN ---- #

    def test_report_false_positive_fixture_passes_end_to_end(self, validator, tmp_path):
        """报告 §八 (vii) 的假阳性夹具：合法输入不得被拦下。

        旧的正则把 ``("A:-3.5":0.25,B:0.25);`` 报成
        ``[CRITICAL] 负分支长度 -3.5``，一份完全合法的输入会被判严重错误。
        """
        tree_file = tmp_path / "quoted_leaf.nwk"
        tree_file.write_text('("A:-3.5":0.25,B:0.25);', encoding="utf-8")

        result = validator.validate_tree_deep(
            tree_file
        )  # 不得抛 NegativeBranchLengthError

        assert result.is_valid is True, result.errors
        assert not any("负分支" in w for w in result.warnings)

    def test_report_false_negative_fixture_is_caught_end_to_end(
        self, validator, tmp_path
    ):
        """报告 §八 (vii) 的假阴性夹具：``:-.5`` 必须真的被拦下（抛退出码 3 的异常）。"""
        from phylodater.core.exceptions import NegativeBranchLengthError

        tree_file = tmp_path / "leading_dot.nwk"
        tree_file.write_text("((A:-.5,B:0.5):0.1,C:0.5);", encoding="utf-8")

        with pytest.raises(NegativeBranchLengthError) as excinfo:
            validator.validate_tree_deep(tree_file)

        assert "负分支长度" in str(excinfo.value)


class TestTaxonomyCircularDependencyValidation:
    """审阅项 C-39：分类学"循环依赖"校验器此前在生产路径零调用点。

    :meth:`DeepValidator.validate_taxonomy_file` 造好了却没人调；这里钉住它
    本身的行为（抛 :class:`TaxonomyConflictError`），并由
    ``tests/unit/cli/test_input_validation.py::TestTaxonomyDeepValidationWired``
    钉住"线已经接到 :class:`InputValidator` 上"。
    """

    @pytest.fixture
    def validator(self):
        return DeepValidator()

    def test_circular_taxonomy_raises_taxonomy_conflict(self, validator, tmp_path):
        from phylodater.core.exceptions import TaxonomyConflictError

        table = tmp_path / "circular.tsv"
        table.write_text(
            "name\ttaxonomy\n" "A\td__Alpha;p__Beta\n" "B\td__Beta;p__Alpha\n",
            encoding="utf-8",
        )

        with pytest.raises(TaxonomyConflictError) as excinfo:
            validator.validate_taxonomy_file(table)

        message = str(excinfo.value)
        assert "循环依赖" in message
        assert "'Alpha' <-> 'Beta'" in message
        # 必须带可执行建议（不能只说"错了"）
        assert excinfo.value.suggestion and "循环引用" in excinfo.value.suggestion

    def test_consistent_taxonomy_table_passes_and_reports_coverage(
        self, validator, tmp_path
    ):
        table = tmp_path / "ok.tsv"
        table.write_text(
            "name\ttaxonomy\n"
            "A\td__Bacteria;p__Firmicutes;c__Bacilli\n"
            "B\td__Archaea;p__Euryarchaeota;c__Methanobacteria\n",
            encoding="utf-8",
        )

        detail = validator.validate_taxonomy_file(table)

        assert detail.is_valid is True
        assert detail.circular_dependencies == []
        assert detail.errors == []
        # 覆盖度必须可核对（B-26 同族：不能"检查了个空"却呈现为通过）
        assert detail.entries_checked == 2

    def test_unparseable_table_discloses_zero_coverage(self, validator, tmp_path):
        """多列表格（没有 d__…;p__… 字符串列）解析不出行 → 必须说明"没查到"。"""
        table = tmp_path / "multi.tsv"
        table.write_text(
            "name\tdomain\tphylum\nhuman\tEukaryota\tChordata\n", encoding="utf-8"
        )

        detail = validator.validate_taxonomy_file(table)

        assert detail.is_valid is True  # 结构没问题，不能因此拦下用户
        assert detail.entries_checked == 0
        assert any("未覆盖任何行" in w for w in detail.warnings), detail.warnings

    def test_missing_and_empty_files_are_errors_not_crashes(self, validator, tmp_path):
        missing = tmp_path / "nope.tsv"
        assert validator.validate_taxonomy_file(missing).is_valid is False

        empty = tmp_path / "empty.tsv"
        empty.write_text("", encoding="utf-8")
        detail = validator.validate_taxonomy_file(empty)
        assert detail.is_valid is False
        assert any("为空" in e for e in detail.errors)


class TestUltrametricityCheck:
    """审阅项 B-18：整个包此前从不检验输入树是否为时间树（超度量）。"""

    @pytest.fixture
    def validator(self):
        return DeepValidator()

    @pytest.fixture
    def temp_dir(self):
        import shutil

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    ULTRAMETRIC = "((A:0.5,B:0.5):0.3,(C:0.4,D:0.4):0.4);"
    NON_ULTRAMETRIC = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"

    def test_ultrametric_tree_is_recognised(self, validator):
        from phylodater.models import PhylogeneticTree

        detail = validator.check_ultrametricity(PhylogeneticTree(self.ULTRAMETRIC))

        assert detail.assessable is True
        assert detail.is_time_tree is True
        assert detail.relative_spread == pytest.approx(0.0)

    def test_non_ultrametric_tree_is_flagged(self, validator):
        from phylodater.models import PhylogeneticTree

        detail = validator.check_ultrametricity(PhylogeneticTree(self.NON_ULTRAMETRIC))

        assert detail.assessable is True
        assert detail.is_time_tree is False
        assert detail.relative_spread > 0.1
        assert "不是" in detail.summarize()
        # 指出具体是哪些叶偏离最远，用户才能定位输入问题
        assert detail.worst_tips and detail.worst_tips[0][0] in {"A", "B", "C", "D"}

    def test_tolerance_is_configurable(self, validator):
        from phylodater.models import PhylogeneticTree

        # root-to-tip: A/B = 0.8, C/D = 0.8001 → 相对极差 ~1.25e-4
        nearly = "((A:0.5,B:0.5):0.3,(C:0.5,D:0.5):0.3001);"

        assert (
            validator.check_ultrametricity(PhylogeneticTree(nearly)).is_time_tree
            is False
        )
        assert (
            validator.check_ultrametricity(
                PhylogeneticTree(nearly), tolerance=1e-3
            ).is_time_tree
            is True
        )

    def test_tree_without_branch_lengths_is_not_assessable(self, validator):
        """纯拓扑树（cladogram）不能声称"是时间树"，也不能冤枉成"不是"。"""
        from phylodater.models import PhylogeneticTree

        detail = validator.check_ultrametricity(PhylogeneticTree("((A,B),(C,D));"))

        assert detail.assessable is False
        assert detail.is_time_tree is False
        assert "无法判定" in detail.summarize()

    def test_accepts_path_and_raw_newick(self, validator, temp_dir):
        tree_file = temp_dir / "t.nwk"
        tree_file.write_text(self.ULTRAMETRIC, encoding="utf-8")

        assert validator.check_ultrametricity(tree_file).is_time_tree is True
        assert validator.check_ultrametricity(self.ULTRAMETRIC).is_time_tree is True
        assert validator.check_ultrametricity(None).assessable is False

    @pytest.mark.parametrize(
        "method,expected",
        [
            ("treepl", "info"),  # 速率平滑：设计上允许非超度量
            ("pathd8", "info"),
            ("r8s", "info"),
            ("lsd2", "warning"),  # 最小二乘定年默认输入已是时间树
            ("iqtree2", "warning"),
            ("mcmctree", "warning"),  # 全局钟把分支长度当时间
            (None, "warning"),
        ],
    )
    def test_severity_is_binned_by_method(self, validator, method, expected):
        assert validator.ultrametric_severity(method) == expected
        # 显式要求时间树时一律升级
        assert validator.ultrametric_severity(method, True) == "error"

    def test_validate_ultrametricity_severity_per_method(self, validator):
        from phylodater.models import PhylogeneticTree

        tree = PhylogeneticTree(self.NON_ULTRAMETRIC)

        errors, warnings, infos = validator.validate_ultrametricity(
            tree, method="treepl"
        )
        assert errors == [] and warnings == [] and len(infos) == 1

        errors, warnings, infos = validator.validate_ultrametricity(tree, method="lsd2")
        assert errors == [] and len(warnings) == 1

        errors, _, _ = validator.validate_ultrametricity(
            tree, method="lsd2", require_time_tree=True
        )
        assert len(errors) == 1

    def test_deep_tree_validation_surfaces_time_tree_warning(self, validator, temp_dir):
        """`--validate` 通路必须真的把这道判定说给用户听（此前零命中）。"""
        tree_file = temp_dir / "nonultra.nwk"
        tree_file.write_text(self.NON_ULTRAMETRIC, encoding="utf-8")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is True  # 只 warn，不拦
        assert any("时间树" in w or "超度量" in w for w in result.warnings)
        assert result.tree_summaries[0]["is_ultrametric"] is False

    def test_ultrametric_tree_summary_marks_true(self, validator, temp_dir):
        tree_file = temp_dir / "ultra.nwk"
        tree_file.write_text(self.ULTRAMETRIC, encoding="utf-8")

        result = validator.validate_tree_deep(tree_file)

        assert result.is_valid is True
        assert result.tree_summaries[0]["is_ultrametric"] is True
        assert not any("不是时间树" in w for w in result.warnings)


class TestIgnoreMalformedOverride:
    """按调用点传入的 ``ignore_malformed`` 必须真正生效。

    ``validate_tree_deep(path, ignore_malformed=...)`` 过去算出了 ``effective_ignore``
    却从不使用，恶意字符检查始终读实例默认值：调用方显式传 ``True`` 仍被硬拒，
    显式传 ``False`` 却被静默降级。ruff 的 F841 报的就是这个死变量。
    """

    @pytest.fixture
    def temp_dir(self):
        import shutil

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def zero_width_tree(self, temp_dir):
        path = temp_dir / "zero_width.nwk"
        # U+200B 零宽空格注入在叶名里
        path.write_text("((A​:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        return path

    def test_per_call_true_downgrades_to_warning(self, zero_width_tree):
        validator = DeepValidator(ignore_malformed=False)  # 实例默认：硬拒

        strict = validator.validate_tree_deep(zero_width_tree)
        overridden = validator.validate_tree_deep(
            zero_width_tree, ignore_malformed=True
        )

        assert strict.is_valid is False
        assert any("零宽" in e for e in strict.errors)
        assert overridden.is_valid is True
        assert not any("零宽" in e for e in overridden.errors)
        assert any("零宽" in w for w in overridden.warnings)

    def test_per_call_false_overrides_permissive_instance(self, zero_width_tree):
        validator = DeepValidator(ignore_malformed=True)  # 实例默认：降级

        result = validator.validate_tree_deep(zero_width_tree, ignore_malformed=False)

        assert result.is_valid is False
        assert any("零宽" in e for e in result.errors)
