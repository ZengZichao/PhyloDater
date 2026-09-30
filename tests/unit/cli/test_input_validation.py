"""输入验证边界条件测试

审阅项 B-24 的行为契约：`validate_calibration_file()` 过去只做五步结构检查
（存在 / 扩展名 / YAML 可解析 / 顶层是 dict / 有 `calibrations` 列表）就返回
`is_valid=True`，CLI 随即打印"校准配置文件格式有效"。本文件里 `TestCalibrationContentValidation`
用报告实测过的五份最普通误写 YAML 逐个钉住新的**内容级**校验。

审阅项 C-39 / C-46 / C-47 / C-48 的契约由本文件后半部分的四个类钉住：
分类学深度校验是否真的接线、未覆盖的格式枚举是否仍默认接受、非 UTF-8 输入是否
真的重试、BioPython 缺失时是否被说成"你的文件有问题"。
"""

import sys
import tempfile
from pathlib import Path
from unittest.mock import Mock

import pytest

from phylodater.cli.main import _clamp_threads, _validate_input_files
from phylodater.services.input_validator import (
    AlignmentFormat,
    InputValidator,
    TreeFormat,
)


class TestInputFileValidation:

    @pytest.fixture
    def mock_logger(self):
        logger = Mock()
        # 让 logger 的方法返回 None（默认行为）
        logger.info.return_value = None
        logger.warning.return_value = None
        logger.error.return_value = None
        logger.success.return_value = None
        return logger

    @pytest.fixture
    def temp_dir(self):
        import shutil

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def valid_files(self, temp_dir):
        tree_file = temp_dir / "tree.nwk"
        tree_file.write_text("((A:0.1,B:0.1):0.2,(C:0.3,D:0.3):0.1);")

        alignment_file = temp_dir / "alignment.fasta"
        alignment_file.write_text(
            ">A\nATCGATCG\n>B\nATCGATCG\n>C\nATCGATCG\n>D\nATCGATCG\n"
        )

        calibrations_file = temp_dir / "calibrations.yaml"
        calibrations_file.write_text(
            "calibrations:\n  - name: Root\n    constraint:\n      type: fixed\n      age: 100\n"
        )

        return tree_file, alignment_file, calibrations_file

    def test_valid_files_pass(self, temp_dir, mock_logger, valid_files):
        tree_file, alignment_file, calibrations_file = valid_files

        parsed_args = Mock()
        parsed_args.tree = tree_file
        parsed_args.sequence = alignment_file
        parsed_args.calibrations = calibrations_file
        parsed_args.taxonomy_file = None

        result = _validate_input_files(parsed_args, mock_logger)

        assert result is True
        mock_logger.error.assert_not_called()

    def test_empty_tree_file_fails(self, temp_dir, mock_logger):
        tree_file = temp_dir / "empty_tree.nwk"
        tree_file.write_text("")

        alignment_file = temp_dir / "alignment.fasta"
        alignment_file.write_text(">A\nATCG\n")

        calibrations_file = temp_dir / "calibrations.yaml"
        calibrations_file.write_text(
            "calibrations:\n  - name: Root\n    constraint:\n      type: fixed\n      age: 100\n"
        )

        parsed_args = Mock()
        parsed_args.tree = tree_file
        parsed_args.sequence = alignment_file
        parsed_args.calibrations = calibrations_file
        parsed_args.taxonomy_file = None

        result = _validate_input_files(parsed_args, mock_logger)

        assert result is False
        mock_logger.error.assert_called()

    def test_empty_alignment_file_fails(self, temp_dir, mock_logger):
        tree_file = temp_dir / "tree.nwk"
        tree_file.write_text("((A:0.1,B:0.1):0.2);")

        alignment_file = temp_dir / "empty_alignment.fasta"
        alignment_file.write_text("")

        calibrations_file = temp_dir / "calibrations.yaml"
        calibrations_file.write_text(
            "calibrations:\n  - name: Root\n    constraint:\n      type: fixed\n      age: 100\n"
        )

        parsed_args = Mock()
        parsed_args.tree = tree_file
        parsed_args.sequence = alignment_file
        parsed_args.calibrations = calibrations_file
        parsed_args.taxonomy_file = None

        result = _validate_input_files(parsed_args, mock_logger)

        assert result is False
        mock_logger.error.assert_called()

    def test_empty_calibration_file_fails(self, temp_dir, mock_logger):
        tree_file = temp_dir / "tree.nwk"
        tree_file.write_text("((A:0.1,B:0.1):0.2);")

        alignment_file = temp_dir / "alignment.fasta"
        alignment_file.write_text(">A\nATCG\n")

        calibrations_file = temp_dir / "empty_calibrations.yaml"
        calibrations_file.write_text("")

        parsed_args = Mock()
        parsed_args.tree = tree_file
        parsed_args.sequence = alignment_file
        parsed_args.calibrations = calibrations_file
        parsed_args.taxonomy_file = None

        result = _validate_input_files(parsed_args, mock_logger)

        assert result is False
        mock_logger.error.assert_called()

    def test_missing_file_fails(self, temp_dir, mock_logger):
        tree_file = temp_dir / "nonexistent_tree.nwk"
        alignment_file = temp_dir / "alignment.fasta"
        alignment_file.write_text(">A\nATCG\n")
        calibrations_file = temp_dir / "calibrations.yaml"
        calibrations_file.write_text(
            "calibrations:\n  - name: Root\n    constraint:\n      type: fixed\n      age: 100\n"
        )

        parsed_args = Mock()
        parsed_args.tree = tree_file
        parsed_args.sequence = alignment_file
        parsed_args.calibrations = calibrations_file
        parsed_args.taxonomy_file = None

        result = _validate_input_files(parsed_args, mock_logger)

        assert result is False
        mock_logger.error.assert_called()

    def test_invalid_tree_format_fails(self, temp_dir, mock_logger):
        tree_file = temp_dir / "invalid_tree.nwk"
        tree_file.write_text("this is not a valid tree")

        alignment_file = temp_dir / "alignment.fasta"
        alignment_file.write_text(">A\nATCG\n")

        calibrations_file = temp_dir / "calibrations.yaml"
        calibrations_file.write_text(
            "calibrations:\n  - name: Root\n    constraint:\n      type: fixed\n      age: 100\n"
        )

        parsed_args = Mock()
        parsed_args.tree = tree_file
        parsed_args.sequence = alignment_file
        parsed_args.calibrations = calibrations_file
        parsed_args.taxonomy_file = None

        result = _validate_input_files(parsed_args, mock_logger)

        assert result is False
        mock_logger.error.assert_called()


class TestLineEndingHandling:
    """换行符归一化（B-24 附带项）。

    注意：``temp_dir`` 只是 :class:`TestInputFileValidation` 的**类内** fixture，
    在这一类里并不可用；pytest 会在 setup 阶段直接报 ``fixture 'temp_dir' not
    found``，使整个测试会话 error。改用内置的 ``tmp_path``。
    """

    def test_crlf_tree_file_converted(self, tmp_path):
        from phylodater.models import PhylogeneticTree

        tree_file = tmp_path / "crlf_tree.nwk"
        tree_file.write_bytes(b"((A:0.1,B:0.1):0.2);\r\n")

        tree = PhylogeneticTree.from_file(tree_file)
        assert tree.num_tips == 2

    def test_lf_file_unchanged(self, tmp_path):
        from phylodater.models import PhylogeneticTree

        tree_file = tmp_path / "lf_tree.nwk"
        tree_file.write_bytes(b"((A:0.1,B:0.1):0.2);\n")

        tree = PhylogeneticTree.from_file(tree_file)
        assert tree.num_tips == 2


class TestCalibrationContentValidation:
    """审阅项 B-24 / B-19：校准 YAML 的**内容级**校验（五份报告实测误写 + 正例）。

    旧实现只做五步结构检查就返回 ``is_valid=True``，CLI 随即打印
    "校准配置文件格式有效: yaml"，而真正_load_时这些点全部被静默丢掉
    （``Loaded 0 calibrations``）。下面把报告里实测过的五种误写逐个钉住：
    校验器必须**拒绝**并给出指向条目序号的错误，CLI 必须返回 False。
    """

    #: 报告 B-24 实测的五种最普通误写（旧实现全部放行）
    MISWRITTEN = {
        "type_at_entry_level": (
            "calibrations:\n"
            "  - name: TestClade\n"
            "    mrca_pair: [human, chimp]\n"
            "    type: uniform\n"
            "    min: 5.0\n"
            "    max: 8.0\n"
        ),
        "min_max_at_entry_level": (
            "calibrations:\n"
            "  - name: TestClade\n"
            "    mrca_pair: [human, chimp]\n"
            "    min: 5.0\n"
            "    max: 8.0\n"
        ),
        "scalar_entry": "calibrations:\n  - TestClade\n",
        "empty_entry": (
            "calibrations:\n"
            "  -\n"
            "  - name: TestClade\n"
            "    constraint:\n"
            "      type: fixed\n"
            "      age: 100.0\n"
        ),
        "missing_constraint": (
            "calibrations:\n" "  - name: TestClade\n" "    mrca_pair: [human, chimp]\n"
        ),
        #: 报告 B-19：``.nan`` 是 YAML 规范字面量，旧实现一路放行到 MCMCTree
        "nan_age": (
            "calibrations:\n"
            "  - name: TestClade\n"
            "    mrca_pair: [human, chimp]\n"
            "    constraint:\n"
            "      type: uniform\n"
            "      min: .nan\n"
            "      max: 100.0\n"
        ),
        "inf_age": (
            "calibrations:\n"
            "  - name: TestClade\n"
            "    mrca_pair: [human, chimp]\n"
            "    constraint:\n"
            "      type: uniform\n"
            "      min: 5.0\n"
            "      max: .inf\n"
        ),
    }

    @pytest.fixture
    def tree_and_alignment(self, tmp_path):
        tree_file = tmp_path / "tree.nwk"
        tree_file.write_text("((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);")
        alignment_file = tmp_path / "alignment.fasta"
        alignment_file.write_text(
            ">human\nATCG\n>chimp\nATCG\n>mouse\nATCG\n>rat\nATCG\n"
        )
        return tree_file, alignment_file

    def _write(self, tmp_path, text, name="calibrations.yaml"):
        path = tmp_path / name
        path.write_text(text, encoding="utf-8")
        return path

    @pytest.mark.parametrize("case", sorted(MISWRITTEN))
    def test_miswritten_yaml_is_rejected_by_validator(self, tmp_path, case):
        result = InputValidator().validate_calibration_file(
            self._write(tmp_path, self.MISWRITTEN[case], f"{case}.yaml")
        )

        assert result.is_valid is False, case
        assert result.errors, case
        # 错误文案必须指向条目序号（报告："报错文案带上条目序号与文件名"）
        assert any("#1" in error or "#2" in error for error in result.errors), case

    @pytest.mark.parametrize("case", sorted(MISWRITTEN))
    def test_cli_does_not_call_miswritten_file_valid(
        self, tmp_path, tree_and_alignment, case
    ):
        tree_file, alignment_file = tree_and_alignment
        parsed_args = Mock()
        parsed_args.tree = tree_file
        parsed_args.sequence = alignment_file
        parsed_args.calibrations = self._write(
            tmp_path, self.MISWRITTEN[case], f"{case}.yaml"
        )
        parsed_args.taxonomy_file = None
        parsed_args.ignore_malformed = False
        parsed_args.skip_length_check = False
        parsed_args.multi_tree_mode = False

        logger = Mock()
        assert _validate_input_files(parsed_args, logger) is False, case
        # 误导文案不得再出现（旧实现在这里打印"校准配置文件格式有效"）
        printed = " ".join(str(call.args[0]) for call in logger.success.call_args_list)
        assert "校准配置文件格式有效" not in printed, case

    def test_wellformed_calibration_file_still_passes(self, tmp_path):
        text = (
            "calibrations:\n"
            "  - name: TestClade\n"
            "    mrca_pair: [human, chimp]\n"
            "    constraint:\n"
            "      type: uniform\n"
            "      min: 5.0\n"
            "      max: 8.0\n"
        )
        result = InputValidator().validate_calibration_file(self._write(tmp_path, text))

        assert result.is_valid is True, result.errors
        assert result.errors == []
        assert any("1" in note for note in result.info)

    def test_empty_calibration_list_is_rejected(self, tmp_path):
        result = InputValidator().validate_calibration_file(
            self._write(tmp_path, "calibrations: []\n")
        )

        assert result.is_valid is False
        assert result.errors

    def test_validator_reports_every_bad_entry_not_just_the_first(self, tmp_path):
        text = (
            "calibrations:\n"
            "  - name: Bad1\n"
            "    constraint:\n"
            "      type: uniform\n"
            "      min: 90.0\n"
            "      max: 10.0\n"
            "  - name: Bad2\n"
            "    constraint:\n"
            "      type: nonsense\n"
            "      age: 5.0\n"
        )
        result = InputValidator().validate_calibration_file(self._write(tmp_path, text))

        assert result.is_valid is False
        assert any("#1" in error for error in result.errors), result.errors
        assert any("#2" in error for error in result.errors), result.errors

    def test_tree_ops_go_through_public_api(self, tmp_path, tree_and_alignment):
        """校验层不得硬依赖 ete3（审阅项 A-5/B-10：Py3.14 上 ete3 无法导入）。"""
        tree_file, alignment_file = tree_and_alignment
        from phylodater.services.deep_validator import DeepValidator

        validator = DeepValidator()
        detail = validator.validate_tree_deep(tree_file)
        assert detail.is_valid
        # B-18：非超度量输入只在 warnings 里出现，不改变 is_valid
        assert "is_ultrametric" in (detail.tree_summaries[0] or {})


class TestClampThreads:
    """测试 --threads/-j 钳制逻辑（#17 杂项修复）"""

    def test_clamps_to_cpu_count_when_exceeded(self):
        """超过物理核数的请求应被钳制到 cpu_count"""
        import os

        max_threads = os.cpu_count() or 1
        assert _clamp_threads(max_threads + 1000) == max_threads

    def test_minimum_is_one(self):
        """下限为 1（负值/0 已在解析阶段被拒绝，此处仅兜底）"""
        assert _clamp_threads(1) == 1
        assert _clamp_threads(0) == 1

    def test_keeps_valid_value(self):
        """区间内的值保持不变"""
        assert _clamp_threads(4) == 4

    def test_accepts_string_int(self):
        """接受字符串形式的整数（argparse 可能传入字符串）"""
        assert _clamp_threads("2") == 2


class TestTaxonomyDeepValidationWired:
    """审阅项 C-39：分类学"循环依赖"校验器必须真的接到生产路径上。

    ``DeepValidator.validate_taxonomy_file``（抛 ``TaxonomyConflictError``、附带
    可执行建议）此前全仓库零调用点：CLI 的 taxonomy 分支只调
    :meth:`InputValidator.validate_taxonomy_file` 的浅检查。现在浅检查通过后会
    再跑一次深度校验，本类钉住"线已接上"的四件事：

    1. 循环表格 → ``is_valid=False``，错误点名"循环依赖"并给出建议；
    2. 验证门的契约是**返回错误清单**而不是抛异常（CLI 不 catch）；
    3. 一致表格 / ``examples/taxonomy.tsv`` 这类合法两列表（首行是表头）仍通过，
       即接线没有引入假阳性；
    4. 多列表格（解析不出分类学字符串）必须披露"本次未覆盖任何行"。
    """

    CIRCULAR = (
        "name\ttaxonomy\n"
        "A\td__Alpha;p__Beta\n"
        "B\td__Beta;p__Alpha\n"  # Beta 既是 Alpha 的子级又是其父级 → 循环
    )
    CONSISTENT = (
        "name\ttaxonomy\n"
        "A\td__Bacteria;p__Firmicutes;c__Bacilli\n"
        "B\td__Archaea;p__Euryarchaeota;c__Methanobacteria\n"
    )

    def _table(self, tmp_path, text, name="taxonomy.tsv"):
        path = tmp_path / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_circular_table_is_rejected_by_validator(self, tmp_path):
        result = InputValidator().validate_taxonomy_file(
            self._table(tmp_path, self.CIRCULAR)
        )

        assert result.is_valid is False
        assert any("循环依赖" in error for error in result.errors), result.errors
        # 可执行建议必须一起来（报告："附带可执行建议"）
        assert any("循环引用" in error for error in result.errors), result.errors

    def test_circular_table_makes_cli_validation_fail(self, tmp_path):
        """CLI 的 taxonomy 分支必须因此返回 False（此前只会打印"格式有效"）。"""
        tree_file = tmp_path / "tree.nwk"
        tree_file.write_text("((A:0.1,B:0.1):0.2,(C:0.3,D:0.3):0.1);", encoding="utf-8")
        alignment_file = tmp_path / "aln.fasta"
        alignment_file.write_text(">A\nATCG\n>B\nATCG\n>C\nATCG\n>D\nATCG\n", "utf-8")
        calibrations_file = tmp_path / "calib.yaml"
        calibrations_file.write_text(
            "calibrations:\n  - name: Root\n"
            "    constraint:\n      type: fixed\n      age: 100\n",
            encoding="utf-8",
        )
        args = Mock()
        args.tree, args.sequence, args.calibrations = (
            tree_file,
            alignment_file,
            calibrations_file,
        )
        args.taxonomy_file = self._table(tmp_path, self.CIRCULAR)
        args.ignore_malformed = False
        args.skip_length_check = False
        args.multi_tree_mode = False
        logger = Mock()

        assert _validate_input_files(args, logger) is False
        printed = " ".join(str(call.args[0]) for call in logger.error.call_args_list)
        assert "循环依赖" in printed
        success = " ".join(str(call.args[0]) for call in logger.success.call_args_list)
        assert "分类学信息文件格式有效" not in success

    def test_validator_actually_delegates_to_deep_validator(
        self, tmp_path, monkeypatch
    ):
        """接线证明：浅检查通过后必须调用一次 DeepValidator 的分类学校验。"""
        from phylodater.services import deep_validator as dv

        calls = []
        original = dv.DeepValidator.validate_taxonomy_file

        def spy(self, file_path):
            calls.append(Path(file_path).name)
            return original(self, file_path)

        monkeypatch.setattr(dv.DeepValidator, "validate_taxonomy_file", spy)
        result = InputValidator().validate_taxonomy_file(
            self._table(tmp_path, self.CONSISTENT)
        )

        assert result.is_valid is True, result.errors
        assert calls == ["taxonomy.tsv"]

    def test_header_two_column_table_is_not_a_false_positive(self, tmp_path):
        """合法的"首行为表头"两列表格必须放行（旧浅检查只看首行第二列）。

        不修这一条，C-39 接上的深度校验对 ``examples/taxonomy.tsv`` 这种
        项目自带格式永远跑不到。
        """
        result = InputValidator().validate_taxonomy_file(
            self._table(tmp_path, self.CONSISTENT)
        )

        assert result.is_valid is True, result.errors
        assert result.format_detected == "taxonomy_table"

    def test_shipped_example_table_still_passes(self):
        """项目自带的 examples/taxonomy.tsv 不能被接线误杀。"""
        example = Path(__file__).resolve().parents[3] / "examples" / "taxonomy.tsv"
        if not example.exists():  # pragma: no cover - 仓库副本缺该文件时跳过
            pytest.skip(f"缺少示例文件 {example}")

        result = InputValidator().validate_taxonomy_file(example)

        assert result.is_valid is True, result.errors

    def test_multi_column_table_discloses_zero_coverage(self, tmp_path):
        """多列表（无 d__… 字符串列）结构合法，但必须说清"这项检查没查到东西"。"""
        result = InputValidator().validate_taxonomy_file(
            self._table(
                tmp_path,
                "name\tdomain\tphylum\nhuman\tEukaryota\tChordata\n",
            )
        )

        assert result.is_valid is True, result.errors
        assert any(
            "未覆盖任何行" in warning for warning in result.warnings
        ), result.warnings

    def test_garbage_two_column_table_is_still_rejected(self, tmp_path):
        """放宽到"表头 + 数据行"不等于接受垃圾文本。"""
        result = InputValidator().validate_taxonomy_file(
            self._table(tmp_path, "a,b\nc,d\ne,f\n")
        )

        assert result.is_valid is False
        assert any("d__" in error for error in result.errors), result.errors


class TestUncoveredFormatIsRejected:
    """审阅项 C-46：格式校验的 if/elif 没有 else → 未覆盖枚举值默认接受。

    ``TreeFormat.STOCKHOLDT`` 与 ``AlignmentFormat.PHYLIP_INTERLEAVED`` 此前没有
    分支，垃圾文本也能拿到 ``is_valid=True``。下面对**每一个**枚举值都跑一份垃圾
    内容，任何一个被默认接受都算失败。
    """

    GARBAGE_TREE = "this is not a tree at all"
    GARBAGE_ALIGNMENT = "this is not an alignment at all"

    def test_garbage_stockholdt_file_is_rejected(self, tmp_path):
        """报告实测夹具：内容全是垃圾文本的 junk.sth。"""
        junk = tmp_path / "junk.sth"
        junk.write_text(self.GARBAGE_TREE, encoding="utf-8")

        result = InputValidator().validate_tree_file(junk)

        assert result.is_valid is False
        assert result.format_detected == "stockholdt"
        assert any("未支持" in error for error in result.errors), result.errors

    @pytest.mark.parametrize("format_", list(TreeFormat))
    def test_every_tree_format_value_is_validated_or_rejected(self, tmp_path, format_):
        """枚举里不能有"悄悄放过"的值——AUTO 也不允许被当成已校验。"""
        junk = tmp_path / "junk.nwk"
        junk.write_text(self.GARBAGE_TREE, encoding="utf-8")

        result = InputValidator()._validate_tree_format(junk, format_)

        assert result.is_valid is False, f"{format_} 被默认接受（无 else 分支）"
        assert result.errors, format_

    @pytest.mark.parametrize("format_", list(AlignmentFormat))
    def test_every_alignment_format_value_is_validated_or_rejected(
        self, tmp_path, format_
    ):
        junk = tmp_path / "junk.fasta"
        junk.write_text(self.GARBAGE_ALIGNMENT, encoding="utf-8")

        result = InputValidator()._validate_alignment_format(junk, format_)

        assert result.is_valid is False, f"{format_} 被默认接受（无 else 分支）"
        assert result.errors, format_

    def test_phylip_interleaved_shares_the_phylip_header_check(self, tmp_path):
        """交错式 PHYLIP 有真的判据（首行 num_taxa seq_length），不是空过。"""
        good = tmp_path / "interleaved.phy"
        good.write_text("2 4\nA ATCG\nB ATCG\n", encoding="utf-8")
        junk = tmp_path / "junk_interleaved.phy"
        junk.write_text(self.GARBAGE_ALIGNMENT, encoding="utf-8")

        validator = InputValidator()
        assert (
            validator._validate_alignment_format(
                good, AlignmentFormat.PHYLIP_INTERLEAVED
            ).is_valid
            is True
        )
        rejected = validator._validate_alignment_format(
            junk, AlignmentFormat.PHYLIP_INTERLEAVED
        )
        assert rejected.is_valid is False
        assert any("PHYLIP" in error for error in rejected.errors), rejected.errors


class TestFormatSignatureTableIsSingleSource:
    """审阅项 C-46 附带项：签名表曾是死数据（FASTA/PHYLIP 从不被读取）。

    判据现在只有 ``*_FORMAT_SIGNATURES`` 一处来源，因此**改表就能改行为**。
    下面逐项篡改表内容并观察检测结果随之变化——这正是"不是死数据"的可证伪证据。
    """

    def test_fasta_detection_follows_the_table(self, tmp_path, monkeypatch):
        aln = tmp_path / "sequences"  # 无扩展名：只能靠内容判据
        aln.write_text(">A\nATCG\n>B\nATCG\n", encoding="utf-8")
        validator = InputValidator()
        assert validator._detect_alignment_format(aln) is AlignmentFormat.FASTA

        monkeypatch.setitem(
            InputValidator.ALIGNMENT_FORMAT_SIGNATURES,
            AlignmentFormat.FASTA,
            [b"~>~"],
        )
        assert validator._detect_alignment_format(aln) is None

    def test_phylip_detection_follows_the_table(self, tmp_path, monkeypatch):
        aln = tmp_path / "alignment"  # 无扩展名
        aln.write_text("2 4\nA ATCG\nB ATCG\n", encoding="utf-8")
        validator = InputValidator()
        assert validator._detect_alignment_format(aln) is AlignmentFormat.PHYLIP

        import re as _re

        monkeypatch.setitem(
            InputValidator.ALIGNMENT_FORMAT_SIGNATURES,
            AlignmentFormat.PHYLIP,
            [_re.compile(rb"^\s*ZZZ\s+\d+\s*$")],
        )
        assert validator._detect_alignment_format(aln) is None
        # 同一个表也驱动**格式验证**（旧实现在 _validate_alignment_format 里又内联
        # 重写了一遍 PHYLIP 正则，于是这里改表对验证毫无影响）
        good = tmp_path / "ok.phy"
        good.write_text("2 4\nA ATCG\nB ATCG\n", encoding="utf-8")
        rejected = validator._validate_alignment_format(good, AlignmentFormat.PHYLIP)
        assert rejected.is_valid is False
        assert any("PHYLIP" in error for error in rejected.errors), rejected.errors

    def test_tree_detection_follows_the_table(self, tmp_path, monkeypatch):
        nexus = tmp_path / "treefile"  # 无扩展名
        nexus.write_text("#NEXUS\nbegin data;\ntree t1 = ((A,B),C);\nend;\n", "utf-8")
        validator = InputValidator()
        assert validator._detect_tree_format(nexus) is TreeFormat.NEXUS

        monkeypatch.setitem(
            InputValidator.TREE_FORMAT_SIGNATURES, TreeFormat.NEXUS, [b"#NOT-A-NEXUS"]
        )
        assert validator._detect_tree_format(nexus) is not TreeFormat.NEXUS

    def test_extension_tables_drive_detection_and_are_not_dead(self, tmp_path):
        """``*_EXTENSIONS`` 此前也是从不被读取的死清单（同 C-46 的形态）。

        现在它们是判据表的派生物，扩展名判定只读 ``*_FORMAT_BY_EXTENSION``——
        包括此前没有任何分支用到的 ``.phylip-interleaved``。
        """
        interleaved = tmp_path / "aln.phylip-interleaved"
        interleaved.write_text("2 4\nA ATCG\nB ATCG\n", encoding="utf-8")

        assert (
            InputValidator()._detect_alignment_format(interleaved)
            is AlignmentFormat.PHYLIP_INTERLEAVED
        )
        assert InputValidator.TREE_EXTENSIONS == frozenset(
            InputValidator.TREE_FORMAT_BY_EXTENSION
        )
        assert InputValidator.ALIGNMENT_EXTENSIONS == frozenset(
            InputValidator.ALIGNMENT_FORMAT_BY_EXTENSION
        )
        # 扩展名 → 格式表与枚举必须一一对应，不能有"表里有、枚举里没有"的值
        assert set(InputValidator.TREE_FORMAT_BY_EXTENSION.values()) <= set(TreeFormat)
        assert set(InputValidator.ALIGNMENT_FORMAT_BY_EXTENSION.values()) <= set(
            AlignmentFormat
        )


class TestEncodingFallbackIsReal:
    """审阅项 C-47：warning 承诺"将尝试其他编码"，但没有任何代码尝试。

    旧行为（报告实测）：一份含 ``0xC9`` 的 latin-1 Newick →
    ``_validate_tree_format`` 给 ``is_valid=True`` + 空话 warning，随后
    ``validate_tree_file`` 整体以原始编解码报错 ``'utf-8' codec can't decode
    byte 0xc9`` 失败。现在按 utf-8-sig → utf-8 → cp1252 → latin-1 **真的**重试，
    并把实际使用的编码记进 info。
    """

    def test_non_utf8_tree_file_is_parsed_after_real_retry(self, tmp_path):
        tree = tmp_path / "latin1.nwk"
        tree.write_bytes("((A\xc9x:0.1,B:0.1):0.2,C:0.3);".encode("latin-1"))

        result = InputValidator().validate_tree_file(tree)

        assert result.is_valid is True, result.errors
        assert result.format_detected == "newick"
        assert not any("codec can't decode" in error for error in result.errors)
        # 报告实际用的编码（既在 warning 也在 info）
        assert any("不是 UTF-8" in w for w in result.warnings), result.warnings
        assert any("cp1252" in i or "latin-1" in i for i in result.info), result.info

    def test_non_utf8_alignment_file_is_parsed_after_real_retry(self, tmp_path):
        aln = tmp_path / "latin1.fasta"
        aln.write_bytes(">A\xc9x\nATCG\n>B\nATCG\n".encode("latin-1"))

        result = InputValidator().validate_alignment_file(aln)

        assert result.is_valid is True, result.errors
        assert any("解码成功" in note for note in result.info), result.info

    def test_best_effort_reader_reports_encoding(self, tmp_path):
        tree = tmp_path / "cp1252.nwk"
        tree.write_bytes("((A:0.1,B:0.1):0.2,C\xe9:0.3);".encode("cp1252"))

        text, encoding = InputValidator()._read_text_best_effort(tree)

        assert encoding in ("cp1252", "latin-1")
        assert "C" in text and "(" in text

    def test_plain_utf8_file_is_not_claimed_to_have_a_bom(self, tmp_path):
        """utf-8-sig 编解码器对无 BOM 文件等价于 utf-8——不得谎报"检测到 BOM"。"""
        tree = tmp_path / "plain.nwk"
        tree.write_text("((A:0.1,B:0.1):0.2,C:0.3);", encoding="utf-8")

        result = InputValidator().validate_tree_file(tree)

        assert result.is_valid is True, result.errors
        assert result.warnings == []
        assert not any("BOM" in note for note in result.info), result.info

    def test_real_bom_is_reported(self, tmp_path):
        tree = tmp_path / "bom.nwk"
        tree.write_bytes(b"\xef\xbb\xbf((A:0.1,B:0.1):0.2,C:0.3);")

        result = InputValidator().validate_tree_file(tree)

        assert result.is_valid is True, result.errors
        assert any("BOM" in note for note in result.info), result.info


class TestBioPythonUnavailableIsBlamedOnEnvironment:
    """审阅项 C-48：BioPython ImportError 被宽 ``except Exception`` 吞掉，
    于是环境装坏时"每个文件都不合法"，且报错让用户去改本来没问题的文件。
    """

    @pytest.fixture
    def without_bio(self, monkeypatch):
        """屏蔽 ``Bio`` 导入（``sys.modules['Bio'] = None`` 让 import 抛 ImportError）。"""
        for name in [n for n in sys.modules if n == "Bio" or n.startswith("Bio.")]:
            monkeypatch.delitem(sys.modules, name, raising=False)
        monkeypatch.setitem(sys.modules, "Bio", None)
        return InputValidator()

    def test_newick_detection_does_not_need_bio(self, tmp_path, without_bio):
        """报告实测夹具：无扩展名的合法 Newick，屏蔽 Bio 后旧实现得到 None。"""
        tree = tmp_path / "no_extension"
        tree.write_text("((A:0.1,B:0.1):0.2,(C:0.3,D:0.3):0.1);", encoding="utf-8")

        assert without_bio._detect_tree_format(tree) is TreeFormat.NEWICK

    def test_tree_parse_error_names_bio(self, tmp_path, without_bio):
        tree = tmp_path / "tree.nwk"
        tree.write_text("((A:0.1,B:0.1):0.2,C:0.3);", encoding="utf-8")

        result = without_bio.validate_tree_file(tree)

        assert result.is_valid is False
        joined = " ".join(result.errors)
        assert "BioPython 不可用" in joined, result.errors
        assert "pip show biopython" in joined, result.errors
        # 不能再把矛头指向用户的文件
        assert "无法自动检测树文件格式" not in joined, result.errors

    def test_alignment_parse_error_names_bio(self, tmp_path, without_bio):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">A\nATCG\n>B\nATCG\n", encoding="utf-8")

        result = without_bio.validate_alignment_file(aln)

        assert result.is_valid is False
        joined = " ".join(result.errors)
        assert "BioPython 不可用" in joined, result.errors
        assert "不是比对文件本身的问题" in joined, result.errors
