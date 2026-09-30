"""
TaxonomyParser 单元测试

- wave-2a: P2 丢级修复（segment 模式下含下划线的分类级值不应被静默丢级）
- B-26: 分类学表整份加载失败时不得伪装成"加载成功 0 条"，一致性检查不得把
  全部叶节点报成"已匹配"
- C-33/C-43/C-44/C-45: 重复等级前缀、重复行、save_to_file(names=…) 数据源、
  一致性清单的确定性排序
"""

from pathlib import Path

import pytest

from phylodater.services.taxonomy_parser import TaxonomyLoadError, TaxonomyParser


class RecordingLogger:
    """记录 PhyloDaterLogger 调用的假日志器（本仓库日志层非 stdlib logging）。"""

    def __init__(self):
        self.records = []

    def _record(self, level, message, **kwargs):
        # 与生产日志层一致的约束：消息必须是字符串（见审阅报告 B-27）
        assert isinstance(
            message, str
        ), f"logger.{level} 收到了非字符串消息 {type(message).__name__}: {message!r}"
        self.records.append((level, message))

    def debug(self, message, **kwargs):
        self._record("DEBUG", message, **kwargs)

    def info(self, message, **kwargs):
        self._record("INFO", message, **kwargs)

    def warning(self, message, **kwargs):
        self._record("WARNING", message, **kwargs)

    def error(self, message, **kwargs):
        self._record("ERROR", message, **kwargs)

    def critical(self, message, **kwargs):
        self._record("CRITICAL", message, **kwargs)

    def logged(self, level, needle):
        return any(lvl == level and needle in msg for lvl, msg in self.records)


def write_table(tmp_path: Path, name: str, text: str) -> Path:
    """写一份分类学表夹具并返回路径。"""
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def instrumented() -> TaxonomyParser:
    """返回一个日志器可被断言的解析器。"""
    parser = TaxonomyParser()
    parser.logger = RecordingLogger()
    return parser


class TestSegmentUnderscoreValues:
    """segment 模式下，值本身含下划线应完整保留，不被误判为新的级别分隔符。"""

    def test_underscore_value_preserved(self):
        """'_d_X_Y_' 中 Y 不是已知级别前缀，整段 'X_Y' 应作为 domain 值保留。"""
        parser = TaxonomyParser(delimiter_mode="segment")
        result = parser.parse("_d_X_Y_")
        assert result == {"domain": "X_Y"}

    def test_underscore_in_multilevel_value(self):
        """phylum 值本身含下划线，应完整保留而非被截断。"""
        parser = TaxonomyParser(delimiter_mode="segment")
        result = parser.parse("GB1_d_Bacteria_p_Proteo_bacteria")
        assert result == {"domain": "Bacteria", "phylum": "Proteo_bacteria"}

    def test_plain_values_still_parsed(self):
        """不含下划线的普通级值仍应正确解析。"""
        parser = TaxonomyParser(delimiter_mode="segment")
        result = parser.parse("_d_Bacteria_p_Cyanobacteriota")
        assert result == {"domain": "Bacteria", "phylum": "Cyanobacteriota"}

    def test_known_prefix_still_delimits(self):
        """已知级别前缀（如 _p_）仍应作为分隔符切分不同级别。"""
        parser = TaxonomyParser(delimiter_mode="segment")
        result = parser.parse("_d_Bacteria_p_Cyanobacteriota_c_")
        assert result.get("domain") == "Bacteria"
        assert result.get("phylum") == "Cyanobacteriota"


class TestLoadFromFileSurfacesFailure:
    """B-26：load_from_file 绝不返回"看起来正常"的空结果。"""

    def test_header_with_domain_but_no_name_column_raises(self, tmp_path):
        """报告中的实测成因：表头有 domain 却没有 name 列。

        旧实现 ``has_header`` 因含 domain 而为真，但 ``row.get("name","")`` 恒为空，
        于是每一行被 ``if not name: continue`` 跳过，只打一条 INFO 就返回 0。
        """
        path = write_table(
            tmp_path,
            "tax.tsv",
            "taxa_id\tdomain\tphylum\n"
            "SpA\tBacteria\tFirmicutes\n"
            "SpB\tArchaea\tEuryarchaeota\n",
        )
        parser = instrumented()

        with pytest.raises(TaxonomyLoadError) as excinfo:
            parser.load_from_file(path)

        # 只捕获 ValueError 的既有调用方（pipeline/cli）同样能接住
        assert isinstance(excinfo.value, ValueError)
        assert "名称列" in str(excinfo.value)
        assert "domain" in str(excinfo.value)  # 报错里列出实际列名
        assert not parser.has_external_file
        assert parser.last_load_report is not None
        assert parser.last_load_report.entries == 0

    def test_no_recognized_rank_columns_raises(self, tmp_path):
        """有名称列但所有分类列名都不被识别 → 报错，而不是静默 0 条。"""
        path = write_table(
            tmp_path,
            "odd.tsv",
            "name\tlin_domain\tlin_phylum\nSpA\tBacteria\tFirmicutes\n",
        )
        parser = instrumented()

        with pytest.raises(TaxonomyLoadError) as excinfo:
            parser.load_from_file(path)
        assert "分类级别列" in str(excinfo.value)

    def test_two_column_strings_that_never_parse_raises(self, tmp_path):
        """两列格式但第二列不是分类串：整份表一行都用不上 → raise。"""
        path = write_table(
            tmp_path,
            "junk.tsv",
            "name\ttaxonomy\nSpA\tBacteria Firmicutes\nSpB\tArchaea Euryarchaeota\n",
        )
        parser = instrumented()

        with pytest.raises(TaxonomyLoadError):
            parser.load_from_file(path)
        report = parser.last_load_report
        assert report.skipped_unparsed == 2
        assert "SpA" in report.failed_names

    def test_header_only_file_raises(self, tmp_path):
        """只有表头、没有数据行 → raise（旧实现返回 0）。"""
        path = write_table(tmp_path, "head.tsv", "name\ttaxonomy\n")
        with pytest.raises(TaxonomyLoadError):
            TaxonomyParser().load_from_file(path)

    def test_empty_file_raises(self, tmp_path):
        """空文件 → raise（TaxonomyLoadError 是 ValueError 子类）。"""
        path = write_table(tmp_path, "empty.tsv", "")
        with pytest.raises(ValueError):
            TaxonomyParser().load_from_file(path)

    def test_blank_lines_only_file_raises(self, tmp_path):
        path = write_table(tmp_path, "blanks.tsv", "\n\n   \n")
        with pytest.raises(TaxonomyLoadError):
            TaxonomyParser().load_from_file(path)

    def test_missing_file_still_filenotfound(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            TaxonomyParser().load_from_file(tmp_path / "nope.tsv")

    def test_failed_load_preserves_previous_table(self, tmp_path):
        """坏文件不得顺手清空已经加载好的表。"""
        good = write_table(
            tmp_path, "good.tsv", "name\ttaxonomy\nSpA\td__Bacteria;p__P\n"
        )
        bad = write_table(tmp_path, "bad.tsv", "taxa_id\tdomain\nX\tBacteria\n")
        parser = instrumented()
        parser.load_from_file(good)

        with pytest.raises(TaxonomyLoadError):
            parser.load_from_file(bad)

        assert parser.has_external_file
        assert sorted(parser._file_loaded_taxonomy) == ["SpA"]
        assert parser.external_file_path == good

    def test_returns_entry_count_and_report(self, tmp_path):
        path = write_table(
            tmp_path,
            "ok.tsv",
            "name\ttaxonomy\nSpA\td__Bacteria;p__P\nSpB\td__Archaea;p__Q\n",
        )
        parser = instrumented()

        assert parser.load_from_file(path) == 2
        assert parser.last_load_report.entries == 2
        assert parser.last_load_report.data_rows == 2
        assert parser.last_load_report.is_success
        assert not parser.last_load_report.has_disclosures

    def test_per_row_skip_reasons_are_tallied(self, tmp_path):
        """B-26：Loaded N 既低估（丢行不报）又高估（空值行照计）→ 现在逐项披露。"""
        path = write_table(
            tmp_path,
            "mixed.tsv",
            "name\ttaxonomy\n"
            "SpA\td__Bacteria;p__P\n"  # 正常
            "SpB\tBacteria;Proteobacteria\n"  # 无前缀 → 解析不出级别
            "SpC\ts__\n"  # GTDB 常见"只有空级"
            "SpA\td__Bacteria;p__OTHER\n"  # 与先入者冲突的重名行
            "SpD;SpE\n",  # 列数不对
        )
        parser = instrumented()

        count = parser.load_from_file(path)

        report = parser.last_load_report
        assert count == 1
        assert report.entries == 1
        assert report.skipped_unparsed == 1
        assert report.skipped_empty_values == 1
        assert report.duplicate_names == 1
        assert report.conflicting_duplicates == 1
        assert report.skipped_malformed == 1
        assert report.data_rows == 5
        assert report.failed_names == ["SpB", "SpC"]
        # 跳过项必须以 WARNING 呈现，而不是只有 INFO
        assert parser.logger.logged("WARNING", "were NOT used")
        assert parser.logger.logged("WARNING", "Duplicate taxonomy entry")

    def test_missing_name_cell_is_tallied_in_multi_column(self, tmp_path):
        """名称列在中间时，空名称行被单列为 missing_name（不当作可用条目）。"""
        path = write_table(
            tmp_path,
            "noname.tsv",
            "domain\tname\tphylum\n"
            "Bacteria\t\tFirmicutes\n"  # 名称列空
            "Archaea\tSpB\tEuryarchaeota\n",
        )
        parser = instrumented()

        assert parser.load_from_file(path) == 1
        report = parser.last_load_report
        assert report.skipped_no_name == 1
        assert "Bacteria" not in parser._file_loaded_taxonomy

    def test_multi_column_by_name_not_position(self, tmp_path):
        """列序无关；表头大小写/前后空格/# 前缀/引号都能识别。"""
        path = write_table(
            tmp_path,
            "cols.tsv",
            '"Phylum", "Name" ,domain\n'
            "Firmicutes,SpA,Bacteria\n"
            "Euryarchaeota,SpB,Archaea\n",
        )
        parser = instrumented()

        assert parser.load_from_file(path) == 2
        assert parser._file_loaded_taxonomy["SpA"] == {
            "domain": "Bacteria",
            "phylum": "Firmicutes",
        }

    def test_gtdb_classification_column_supported(self, tmp_path):
        """GTDB metadata 风格：名称列 + 整条 classification 串。"""
        path = write_table(
            tmp_path,
            "gtdb.tsv",
            "#taxon\tclassification\tassembly_level\n"
            "g__Escherichia;s__\td__Bacteria;p__Proteobacteria;g__Escherichia;s__\tWGS\n",
        )
        parser = instrumented()

        assert parser.load_from_file(path) == 1
        assert (
            parser._file_loaded_taxonomy["g__Escherichia;s__"]["domain"] == "Bacteria"
        )

    def test_extra_rank_columns_are_not_dropped(self, tmp_path):
        """strain/kingdom 列同样是分类信息，旧实现按固定 7 列读会静默丢掉。"""
        path = write_table(
            tmp_path,
            "ranks.tsv",
            "name\tdomain\tkingdom\tstrain\nSpA\tBacteria\tnbac\tH1\n",
        )
        parser = instrumented()

        parser.load_from_file(path)
        assert parser._file_loaded_taxonomy["SpA"] == {
            "domain": "Bacteria",
            "kingdom": "nbac",
            "strain": "H1",
        }

    def test_trailing_empty_cell_does_not_drop_row(self, tmp_path):
        """尾部空列在真实 TSV 里极常见：按空值处理，不能把整行丢掉。"""
        path = write_table(
            tmp_path,
            "ragged.tsv",
            "name\tdomain\tphylum\tstrain\n"
            "SpA\tBacteria\tFirmicutes\tH1\n"
            "SpB\tArchaea\tEuryarchaeota\t\n",
        )
        parser = instrumented()

        assert parser.load_from_file(path) == 2
        assert parser.last_load_report.short_rows == 1
        assert parser.logger.logged("WARNING", "fewer columns than the header")

    def test_bom_and_tab_still_supported(self, tmp_path):
        """回归：BOM 表头 + 两列 TSV 仍然按名称对应。"""
        path = tmp_path / "bom.tsv"
        path.write_bytes("﻿name\ttaxonomy\nSpA\td__Bacteria;p__P\n".encode("utf-8"))
        parser = instrumented()

        assert parser.load_from_file(path) == 1
        assert "SpA" in parser._file_loaded_taxonomy

    def test_illegal_control_character_still_wrapped_as_valueerror(self, tmp_path):
        """§六.56 确认过的行为：解析异常包成带原因的 ValueError 上抛。"""
        path = tmp_path / "ctrl.tsv"
        path.write_text("name\ttaxonomy\nSpA\td__Bac\x01teria\n", encoding="utf-8")
        parser = instrumented()

        with pytest.raises(ValueError) as excinfo:
            parser.load_from_file(path)
        assert not isinstance(excinfo.value, TaxonomyLoadError)


class TestCheckLabelConsistencyIsNotDeceptive:
    """B-26：没有可用分类表时，一致性检查绝不能报"全部匹配"。"""

    def test_returns_none_when_no_table_loaded(self):
        parser = instrumented()

        result = parser.check_label_consistency(["SpA", "SpB", "SpC"])

        assert result is None
        assert parser.logger.logged("WARNING", "NOT PERFORMED")

    def test_report_marks_check_as_not_performed(self, tmp_path):
        parser = instrumented()

        report = parser.evaluate_label_consistency(["SpA", "SpB"])

        assert report.performed is False
        assert report.matched == []
        assert report.unmatched_in_tree == []
        assert not report.fully_matched
        assert "未执行" in report.summary()

    def test_partial_match_reports_both_unmatched_lists(self, tmp_path):
        path = write_table(
            tmp_path, "t.tsv", "name\ttaxonomy\nSpA\td__Bacteria\nSpD\td__Archaea\n"
        )
        parser = instrumented()
        parser.load_from_file(path)

        matched, unmatched_tree, unmatched_table = parser.check_label_consistency(
            ["SpC", "SpB", "SpA"]
        )

        assert matched == ["SpA"]
        assert unmatched_tree == ["SpB", "SpC"]
        assert unmatched_table == ["SpD"]

    def test_zero_overlap_is_reported_not_fabricated(self, tmp_path):
        path = write_table(tmp_path, "t.tsv", "name\ttaxonomy\nOther\td__Bacteria\n")
        parser = instrumented()
        parser.load_from_file(path)

        matched, unmatched_tree, unmatched_table = parser.check_label_consistency(
            ["SpA", "SpB"]
        )

        assert matched == []
        assert unmatched_tree == ["SpA", "SpB"]
        assert unmatched_table == ["Other"]
        assert parser.logger.logged("WARNING", "0 of 2 tree tip(s)")

    def test_lists_are_sorted_for_reproducibility(self, tmp_path):
        """C-45：list(set) 顺序随 PYTHONHASHSEED 变化 → 三个清单必须已排序。"""
        names = [f"Taxon_{i:03d}" for i in range(200)]
        path = write_table(
            tmp_path,
            "big.tsv",
            "name\ttaxonomy\n" + "".join(f"{n}\td__Bacteria\n" for n in names),
        )
        parser = instrumented()
        parser.load_from_file(path)

        matched, unmatched_tree, unmatched_table = parser.check_label_consistency(
            names + ["Zz_A", "Aa_B"]
        )

        assert matched == sorted(names)
        assert unmatched_tree == ["Aa_B", "Zz_A"]
        assert unmatched_table == []
        # 两次调用给出完全相同的顺序
        assert parser.check_label_consistency(names + ["Zz_A", "Aa_B"]) == (
            matched,
            unmatched_tree,
            unmatched_table,
        )

    def test_consistency_after_failed_load_is_none(self, tmp_path):
        """表加载失败 + 一致性检查：两步都得说"没做"，而不是 100% 通过。"""
        bad = write_table(tmp_path, "bad.tsv", "taxa_id\tdomain\nX\tBacteria\n")
        parser = instrumented()
        with pytest.raises(TaxonomyLoadError):
            parser.load_from_file(bad)

        assert parser.check_label_consistency(["SpA", "SpB", "SpC"]) is None


class TestDuplicateRowPolicy:
    """C-43：同名重复行不再"后写者静默胜出"。"""

    def test_first_row_wins_and_conflict_is_warned(self, tmp_path):
        path = write_table(
            tmp_path,
            "dup.tsv",
            "name\ttaxonomy\nSpA\td__Bacteria;p__P1\nSpA\td__Bacteria;p__P2\n",
        )
        parser = instrumented()

        count = parser.load_from_file(path)

        assert count == 1
        assert parser._file_loaded_taxonomy["SpA"]["phylum"] == "P1"
        assert parser.last_load_report.conflicting_duplicates == 1
        assert parser.logger.logged("WARNING", "Duplicate taxonomy entry")

    def test_identical_duplicate_is_merged_quietly(self, tmp_path):
        path = write_table(
            tmp_path,
            "dup2.tsv",
            "name\ttaxonomy\nSpA\td__Bacteria;p__P1\nSpA\td__Bacteria;p__P1\n",
        )
        parser = instrumented()

        assert parser.load_from_file(path) == 1
        assert parser.last_load_report.duplicate_names == 1
        assert parser.last_load_report.conflicting_duplicates == 0
        assert not parser.logger.logged("WARNING", "Duplicate taxonomy entry")

    def test_duplicate_count_is_not_reported_as_entries(self, tmp_path):
        """旧实现 loaded_count 按行计，3 行 2 条会被报成 3。"""
        path = write_table(
            tmp_path,
            "dup3.tsv",
            "name\ttaxonomy\n"
            "SpA\td__Bacteria;p__P1\n"
            "SpA\td__Bacteria;p__P2\n"
            "SpB\td__Archaea;p__P3\n",
        )
        parser = instrumented()

        count = parser.load_from_file(path)

        assert count == 2 == len(parser._file_loaded_taxonomy)
        assert parser.last_load_report.data_rows == 3


class TestSaveToFileNameSource:
    """C-44：save_to_file(names=…) 必须看得见刚加载的分类表。"""

    def test_names_argument_finds_table_entries(self, tmp_path):
        path = write_table(
            tmp_path,
            "t.tsv",
            "name\ttaxonomy\nSpA\td__Bacteria;p__P\nSpB\td__Archaea;p__Q\n",
        )
        parser = instrumented()
        parser.load_from_file(path)
        out = tmp_path / "out.tsv"

        written = parser.save_to_file(out, names=["SpA", "SpB"])

        assert written == 2
        assert out.exists()
        lines = out.read_text(encoding="utf-8").strip().splitlines()
        assert len(lines) == 3
        assert lines[1].startswith("SpA\t")

    def test_missing_names_are_disclosed_and_counted(self, tmp_path):
        path = write_table(tmp_path, "t.tsv", "name\ttaxonomy\nSpA\td__Bacteria\n")
        parser = instrumented()
        parser.load_from_file(path)
        out = tmp_path / "out.tsv"

        written = parser.save_to_file(out, names=["SpA", "Ghost1", "Ghost2"])

        assert written == 1
        assert parser.logger.logged("WARNING", "2 of 3 requested name(s)")

    def test_nothing_found_returns_zero_and_writes_no_file(self, tmp_path):
        parser = instrumented()
        out = tmp_path / "out.tsv"

        assert parser.save_to_file(out, names=["Ghost"]) == 0
        assert not out.exists()
        assert parser.logger.logged("WARNING", "was NOT written")

    def test_cache_only_names_still_saved(self, tmp_path):
        parser = instrumented()
        assert parser.parse("d__Bacteria;p__Proteobacteria") is not None
        out = tmp_path / "out.tsv"

        assert parser.save_to_file(out) == 1
        assert out.exists()


class TestDuplicateRankPrefixDisclosure:
    """C-33：同一等级前缀在名称里出现两次时，覆盖必须留痕。"""

    @pytest.mark.parametrize("mode", ["greedy", "reverse", "segment"])
    def test_repeated_rank_prefix_warns(self, mode):
        parser = TaxonomyParser(delimiter_mode=mode)
        parser.logger = RecordingLogger()

        result = parser.parse("GB1_d_Bacteria_g_GenusA_g_GenusB")

        assert result["genus"] == "GenusB"  # 语义不变：后值胜出
        assert parser.logger.logged("WARNING", "同一等级前缀重复出现")
        assert parser.logger.logged("WARNING", "GenusA")

    def test_no_warning_without_repetition(self):
        parser = TaxonomyParser()
        parser.logger = RecordingLogger()

        assert parser.parse("GB1_d_Bacteria_p_Cyanobacteriota_g_GenusX")
        assert not any(level == "WARNING" for level, _ in parser.logger.records)

    def test_taxonomy_string_with_repeated_prefix_warns(self):
        parser = TaxonomyParser()
        parser.logger = RecordingLogger()

        result = parser._parse_taxonomy_string("d__X;d__Y")

        assert result == {"domain": "Y"}
        assert parser.logger.logged("WARNING", "同一等级前缀重复出现")

    def test_single_underline_string_with_repeated_prefix_warns(self):
        parser = TaxonomyParser()
        parser.logger = RecordingLogger()

        result = parser._parse_taxonomy_string("d_X;d_Y")

        assert result == {"domain": "Y"}
        assert parser.logger.logged("WARNING", "同一等级前缀重复出现")

    def test_reverse_mode_agreement_warning_preserved(self):
        """§六.45 的正确范式不能被回归：两方向一致/不一致都要有 warning。"""
        parser = TaxonomyParser(delimiter_mode="reverse")
        parser.logger = RecordingLogger()

        result = parser.parse("GB1_d_Bacteria_s_sp")

        assert result == {"domain": "Bacteria", "species": "sp"}
        assert parser.logger.logged("WARNING", "Ambiguous prefix detected")


class TestTableReloadInvalidatesParseCache:
    """表格来源的解析结果被缓存后，换表必须让缓存失效。"""

    def test_second_load_wins(self, tmp_path):
        first = write_table(
            tmp_path, "first.tsv", "name\ttaxonomy\nSpA\td__Bacteria;p__P1\n"
        )
        second = write_table(
            tmp_path, "second.tsv", "name\ttaxonomy\nSpA\td__Archaea;p__P2\n"
        )
        parser = instrumented()

        parser.load_from_file(first)
        assert parser.parse("SpA")["phylum"] == "P1"
        parser.load_from_file(second)
        assert parser.parse("SpA")["phylum"] == "P2"
