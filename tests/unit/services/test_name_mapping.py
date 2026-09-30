"""
NameMapping 服务单元测试

测试名称映射服务的各项功能
"""

import shutil
import tempfile
from pathlib import Path

import pytest

from phylodater.services.name_mapping import (
    NameCollisionError,
    NameMappingEntry,
    NameMappingManager,
    NameShortenMode,
    resolve_short_name,
    sanitize_names_uniquely,
)


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


class TestNameMappingManager:
    """NameMappingManager 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="name_mapping_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def mapper(self):
        """创建映射器实例"""
        return NameMappingManager()

    def test_add_entry(self, mapper):
        """测试添加映射条目"""
        short_name = mapper.add_entry("original_name", NameShortenMode.FULL)

        assert short_name == "original_name"
        assert "original_name" in mapper._entries

    def test_get_short_name(self, mapper):
        """测试获取简化名称"""
        mapper.add_entry("seq1", NameShortenMode.FULL)

        short = mapper.get_short_name("seq1")

        assert short is not None

    def test_get_short_name_not_found(self, mapper):
        """测试获取不存在的简化名称"""
        short = mapper.get_short_name("nonexistent")

        assert short is None

    def test_get_original_name(self, mapper):
        """测试获取原始名称"""
        mapper.add_entry("seq1", NameShortenMode.FULL)

        original = mapper.get_original_name("seq1")

        assert original == "seq1"

    def test_get_original_name_not_found(self, mapper):
        """测试获取不存在的原始名称"""
        original = mapper.get_original_name("nonexistent")

        assert original is None

    def test_get_mapping_dict(self, mapper):
        """测试获取映射字典"""
        mapper.add_entry("seq1", NameShortenMode.FULL)
        mapper.add_entry("seq2", NameShortenMode.FULL)

        mapping = mapper.get_mapping_dict()

        assert len(mapping) == 2
        assert "seq1" in mapping

    def test_get_reverse_mapping_dict(self, mapper):
        """测试获取反向映射字典"""
        mapper.add_entry("seq1", NameShortenMode.FULL)
        mapper.add_entry("seq2", NameShortenMode.FULL)

        reverse = mapper.get_reverse_mapping_dict()

        assert len(reverse) == 2

    def test_normalize_name(self, mapper):
        """测试名称规范化"""
        # NameMappingManager有sanitize_name静态方法
        assert NameMappingManager.sanitize_name("Test Name") is not None

    def test_serial_mode(self, mapper):
        """测试串行编号模式"""
        short1 = mapper.add_entry("seq1", NameShortenMode.SERIAL)
        short2 = mapper.add_entry("seq2", NameShortenMode.SERIAL)

        assert short1 != short2
        assert len(short1) <= 10

    def test_accession_mode(self, mapper):
        """测试登录号模式"""
        short = mapper.add_entry("ABCD1234.1 Homo sapiens", NameShortenMode.ACCESSION)

        assert short is not None
        assert len(short) <= 30

    def test_save_and_load_file(self, mapper, temp_dir):
        """测试保存和加载文件"""
        mapper.add_entry("seq1", NameShortenMode.FULL)
        mapper.add_entry("seq2", NameShortenMode.FULL)

        file_path = temp_dir / "mapping.txt"
        mapper.save_to_file(file_path)

        assert file_path.exists()

        # 创建新的mapper并加载
        new_mapper = NameMappingManager()
        new_mapper.load_from_file(file_path)

        assert new_mapper.get_short_name("seq1") is not None


class TestNameMappingEntry:
    """NameMappingEntry 测试类"""

    def test_create_entry(self):
        """测试创建映射条目"""
        entry = NameMappingEntry(
            original="Homo sapiens",
            short="homo_sap",
            accession="ABCD1234",
            taxonomy={"genus": "Homo", "species": "sapiens"},
            method=NameShortenMode.ACCESSION,
        )

        assert entry.original == "Homo sapiens"
        assert entry.short == "homo_sap"
        assert entry.accession == "ABCD1234"


class TestNameMappingEdgeCases:
    """NameMapping 边界情况测试"""

    @pytest.fixture
    def mapper(self):
        return NameMappingManager()

    def test_empty_name(self, mapper):
        """测试空名称"""
        with pytest.raises(ValueError):
            mapper.add_entry("", NameShortenMode.FULL)

    def test_none_name(self, mapper):
        """测试None名称"""
        with pytest.raises(TypeError):
            mapper.add_entry(None, NameShortenMode.FULL)

    def test_duplicate_mapping(self, mapper):
        """测试重复映射"""
        short1 = mapper.add_entry("seq1", NameShortenMode.FULL)
        short2 = mapper.add_entry("seq1", NameShortenMode.FULL)  # 返回已存在的

        assert short1 == short2

    def test_many_entries(self, mapper):
        """测试多个条目"""
        for i in range(100):
            mapper.add_entry(f"seq{i}", NameShortenMode.SERIAL)

        assert len(mapper._entries) == 100

    def test_extract_accession(self, mapper):
        """测试提取登录号"""
        accession = mapper._extract_accession("ABCD1234.1 Homo sapiens")

        assert accession is not None

    def test_parse_taxonomy(self, mapper):
        """测试解析分类学信息"""
        taxonomy = mapper._parse_taxonomy("Homo sapiens")

        assert taxonomy is not None


class TestNameMappingCache:
    """测试实例级解析缓存（#16 重复代码收敛 - name_mapping）"""

    @pytest.fixture
    def mapper(self):
        return NameMappingManager()

    def test_extract_accession_is_cached(self, mapper):
        """重复调用 _extract_accession 应命中缓存，返回一致结果"""
        name = "GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota"
        first = mapper._extract_accession(name)
        # 缓存应已写入
        assert name in mapper._accession_cache
        assert mapper._accession_cache[name] == first

        second = mapper._extract_accession(name)
        assert second == first
        # 缓存条目数不应因重复调用而增长
        assert len(mapper._accession_cache) == 1

    def test_parse_taxonomy_is_cached(self, mapper):
        """重复调用 _parse_taxonomy 应命中缓存，返回一致结果"""
        name = "_d_Bacteria_p_Cyanobacteriota_c_Algae_o_Unknown"
        first = mapper._parse_taxonomy(name)
        assert name in mapper._taxonomy_cache
        assert mapper._taxonomy_cache[name] == first

        second = mapper._parse_taxonomy(name)
        assert second == first
        assert len(mapper._taxonomy_cache) == 1

    def test_add_entry_reuses_parsing_cache(self, mapper):
        """相同名称重复 add_entry 不应重复解析（缓存命中）"""
        name = "GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota"
        mapper.add_entry(name, NameShortenMode.ACCESSION)
        mapper.add_entry(name, NameShortenMode.ACCESSION)  # 返回已存在条目

        # 解析缓存仍只保留一条
        assert len(mapper._accession_cache) == 1
        assert len(mapper._taxonomy_cache) == 1


class TestShortNameUniqueness:
    """B-21：共享登录号的两个序列必须得到两个不同短名，且反查不丢分类单元。

    这是本仓库唯一一条会静默改变"输出树里有谁"的缺陷路径：MCMCTree 通路用
    登录号当叶短名（``mcmctree_method._build_name_mapping`` →
    ``add_entry(method=ACCESSION)``），随后 ``_restore_leaf_names`` 用反查表把
    短名换回全名。反查表一旦遮蔽，一个分类单元从输出树消失、另一个出现两次，
    而输出树看上去完全正常。
    """

    # 同一 assembly（同一 GCA_…）拆成两条分片：GTDB 风格数据里很常见
    SHARED_ACCESSION_A = (
        "GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_c_Algo_f_X_g_Y_s_u"
    )
    SHARED_ACCESSION_B = (
        "GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_c_Algo_f_X_g_Y_s_different"
    )

    @pytest.fixture
    def mapper(self):
        manager = NameMappingManager()
        manager.logger = RecordingLogger()
        return manager

    def _build(self, mapper):
        first = mapper.add_entry(self.SHARED_ACCESSION_A, NameShortenMode.ACCESSION)
        second = mapper.add_entry(self.SHARED_ACCESSION_B, NameShortenMode.ACCESSION)
        return first, second

    def test_shared_accession_yields_distinct_short_names(self, mapper):
        """核心复现：两条共享登录号的序列得到不同短名。"""
        first, second = self._build(mapper)

        assert first != second
        assert first == "GB_GCA_000252485_1"
        assert second == "GB_GCA_000252485_1_2"
        # 登录号本身仍然如实记录（两条都指向同一 assembly）
        assert mapper.get_accession(self.SHARED_ACCESSION_A) == "GB_GCA_000252485_1"
        assert mapper.get_accession(self.SHARED_ACCESSION_B) == "GB_GCA_000252485_1"

    def test_reverse_mapping_keeps_both_taxa(self, mapper):
        """反查表必须两个 taxon 都在，且各自还原回自己（不是同一个全名）。"""
        first, second = self._build(mapper)
        reverse = mapper.get_reverse_mapping_dict()

        assert len(reverse) == 2
        assert reverse[first] == self.SHARED_ACCESSION_A
        assert reverse[second] == self.SHARED_ACCESSION_B
        # 复现 _restore_leaf_names 的用法：dict 推导反查不丢任何一行
        forward = mapper.get_mapping_dict()
        assert len(set(forward.values())) == len(forward)
        assert {short: full for full, short in forward.items()} == reverse

    def test_collision_is_disclosed_not_silent(self, mapper):
        """消歧必须留下痕迹（WARNING），不能静默改名。"""
        self._build(mapper)

        assert mapper.logger.logged("WARNING", "短名冲突已消歧")
        assert mapper.logger.logged("WARNING", "GB_GCA_000252485_1_2")

    def test_many_sequences_sharing_one_accession(self, mapper):
        """N 条共享登录号 → N 个互不相同的短名，一个都不能少。"""
        names = [f"GB_GCA_000999999.1_d_Bacteria_p_P_g_G_s_part{i}" for i in range(5)]
        shorts = [mapper.add_entry(n, NameShortenMode.ACCESSION) for n in names]

        assert len(set(shorts)) == 5
        assert len(mapper.get_reverse_mapping_dict()) == 5
        for name, short in zip(names, shorts):
            assert mapper.get_original_name(short) == name

    def test_disambiguated_names_respect_length_limits(self, mapper):
        """消歧后缀必须落在本模式的长度预算内（SERIAL≤10 / ACCESSION≤30）。"""
        long_accession_names = [
            f"GB_GCA_{'1' * 24}.{i}_d_Bacteria" for i in range(1, 4)
        ]
        shorts = [
            mapper.add_entry(n, NameShortenMode.ACCESSION) for n in long_accession_names
        ]
        assert len(set(shorts)) == 3
        assert all(len(s) <= 30 for s in shorts)

        serial_mapper = NameMappingManager()
        serial_mapper.logger = RecordingLogger()
        # SERIAL 短名靠计数器本就唯一；此处校验"计数器与文件里既有短名撞车"时仍唯一
        serial_mapper._short_to_original["S00001"] = "someone-else"
        got = [
            serial_mapper.add_entry(n, NameShortenMode.SERIAL)
            for n in ("a1", "a2", "a3")
        ]
        assert len(set(got)) == 3
        assert all(len(s) <= 10 for s in got)
        assert "S00001" not in got

    def test_restore_leaf_names_flow_preserves_every_tip(self, mapper):
        """端到端复现 MCMCTree 通路：全名→短名→（树被换成短名）→还原全名。

        旧实现把两个叶都还原成同一个全名：一个分类单元从输出树消失、另一个出现
        两次，而输出树看上去完全正常（这是本仓库唯一一条静默改变"树里有谁"的缺陷）。
        """
        tips = [self.SHARED_ACCESSION_A, self.SHARED_ACCESSION_B]
        for tip in tips:
            mapper.add_entry(tip, NameShortenMode.ACCESSION)

        # _build_name_mapping 的产物（mcmctree_method / treepl_method）
        forward = mapper.get_mapping_dict()
        assert len(set(forward.values())) == len(forward)

        # 输入侧：树与比对都换成短名
        shortened = [forward[tip] for tip in tips]
        assert len(set(shortened)) == 2

        # 输出侧：_restore_leaf_names 用反向表还原
        reverse = {short: full for full, short in forward.items()}
        restored = [reverse[short] for short in shortened]
        assert sorted(restored) == sorted(tips)
        assert len(set(restored)) == 2

    def test_integrity_helpers_report_no_problem_when_unique(self, mapper):
        self._build(mapper)

        assert mapper.find_duplicate_short_names() == {}
        assert mapper.check_mapping_integrity() == []
        mapper.assert_mapping_integrity()  # 不抛异常

    def test_reverse_lookup_refuses_silently_shadowed_entry(self, mapper):
        """反查表若会丢掉某个 taxon，默认必须 raise 而不是给出"看起来正常"的表。"""
        self._build(mapper)
        # 模拟外部/旧文件造成的破坏：两条原始名共用一个短名
        entry = NameMappingEntry(
            original="ghost_taxon",
            short="GB_GCA_000252485_1",
            accession="GB_GCA_000252485_1",
            taxonomy={},
            method=NameShortenMode.ACCESSION,
        )
        mapper._entries["ghost_taxon"] = entry

        with pytest.raises(NameCollisionError):
            mapper.get_reverse_mapping_dict()
        assert mapper.check_mapping_integrity()
        with pytest.raises(NameCollisionError):
            mapper.assert_mapping_integrity()

        # strict=False 只降级为 warning，用于诊断场景
        recovered = mapper.get_reverse_mapping_dict(strict=False)
        assert isinstance(recovered, dict)
        assert mapper.logger.logged("WARNING", "名称映射完整性检查失败")

    def test_load_from_file_refuses_duplicate_short_names(self, mapper, tmp_path):
        """映射文件里两行共用同一短名 → 拒绝加载（B-21 的持久化形态）。"""
        bad = tmp_path / "mapping_bad.tsv"
        bad.write_text(
            "original\tshort\taccession\ttaxonomy\tmethod\n"
            "taxon_A\tS1\tACC_A\t\tfull\n"
            "taxon_B\tS1\tACC_B\t\tfull\n",
            encoding="utf-8",
        )

        good = tmp_path / "mapping_good.tsv"
        good.write_text(
            "original\tshort\taccession\ttaxonomy\tmethod\n"
            "taxon_A\tS1\tACC_A\tdomain=Bacteria\tfull\n",
            encoding="utf-8",
        )
        mapper.load_from_file(good)

        with pytest.raises(NameCollisionError):
            mapper.load_from_file(bad)
        # 加载坏文件不得清空已有映射
        assert mapper.get_short_name("taxon_A") == "S1"

    def test_load_from_file_refuses_inconsistent_original_rows(self, mapper, tmp_path):
        """同一 original 出现两次且短名不同 → 还原目标不唯一，直接拒绝。"""
        bad = tmp_path / "mapping_dup_original.tsv"
        bad.write_text(
            "original\tshort\taccession\ttaxonomy\tmethod\n"
            "taxon_A\tS1\t\t\tfull\n"
            "taxon_A\tS2\t\t\tfull\n",
            encoding="utf-8",
        )

        with pytest.raises(NameCollisionError):
            mapper.load_from_file(bad)

    def test_saved_file_round_trips_disambiguated_names(self, mapper, tmp_path):
        """消歧后的短名可原样写出并读回，计数器不会退回已用编号。"""
        first, second = self._build(mapper)
        path = tmp_path / "mapping.tsv"
        mapper.save_to_file(path)

        reloaded = NameMappingManager()
        reloaded.logger = RecordingLogger()
        reloaded.load_from_file(path)

        assert reloaded.get_mapping_dict() == {
            self.SHARED_ACCESSION_A: first,
            self.SHARED_ACCESSION_B: second,
        }
        assert reloaded.get_reverse_mapping_dict() == mapper.get_reverse_mapping_dict()

    def test_accession_lookup_returns_all_sharers(self, mapper):
        """按登录号反查必须是多值表，不能只显示"最后一个"或"第一个"。"""
        self._build(mapper)

        assert mapper.get_originals_by_accession("GB_GCA_000252485_1") == [
            self.SHARED_ACCESSION_A,
            self.SHARED_ACCESSION_B,
        ]
        assert mapper.get_accession(self.SHARED_ACCESSION_B) == "GB_GCA_000252485_1"


class TestResolveShortName:
    """``resolve_short_name``：短名唯一性的单一实现点。"""

    def test_free_slot_is_returned_unchanged(self):
        assert resolve_short_name("ABC", {}, "owner") == "ABC"

    def test_same_owner_is_idempotent(self):
        claimed = {"ABC": "owner"}
        assert resolve_short_name("ABC", claimed, "owner") == "ABC"

    def test_other_owner_gets_counter_suffix(self):
        claimed = {"ABC": "other"}
        assert resolve_short_name("ABC", claimed, "owner") == "ABC_2"

    def test_suffixes_within_length_budget(self):
        claimed = {"ABCDEFGH": "a", "ABCDEFGH_2": "b"}
        got = resolve_short_name("ABCDEFGH", claimed, "c", max_length=10)
        assert got == "ABCDEFGH_3"
        assert len(got) <= 10

    def test_truncation_happens_before_suffixing(self):
        got = resolve_short_name("A" * 30, {"A" * 30: "other"}, "owner", max_length=30)
        assert got == "A" * 28 + "_2"
        assert len(got) == 30

    def test_impossible_budget_raises_instead_of_looping(self):
        with pytest.raises(NameCollisionError):
            resolve_short_name("ab", {"a": "x", "_": "y"}, "owner", max_length=1)

    def test_empty_base_raises(self):
        with pytest.raises(NameCollisionError):
            resolve_short_name("", {}, "owner")


class TestSanitizeNameUniqueness:
    """C-30：``sanitize_name`` 截断会让不同全名塌成同一短名，必须有消歧入口。"""

    COLLIDING = [
        "Escherichia_coli_strain_ABCDEFGHIJ_more_sequence_data",
        "Escherichia_coli_strain_ABCDEFGHIJ_other_sequence_data",
    ]

    def test_plain_sanitize_still_collides(self):
        """纯函数语义保持不变（单名清洗、不消歧）。"""
        first, second = self.COLLIDING
        assert NameMappingManager.sanitize_name(first) == (
            NameMappingManager.sanitize_name(second)
        )

    def test_batch_entry_disambiguates(self):
        mapping = sanitize_names_uniquely(self.COLLIDING)

        assert len(set(mapping.values())) == 2
        assert all(len(s) <= 30 for s in mapping.values())
        assert mapping[self.COLLIDING[0]] == "Escherichia_coli_strain_ABCDEF"
        assert mapping[self.COLLIDING[1]] != mapping[self.COLLIDING[0]]

    def test_batch_entry_is_deterministic(self):
        """同一输入顺序 → 同一输出（可复现性要求）。"""
        assert sanitize_names_uniquely(self.COLLIDING) == sanitize_names_uniquely(
            self.COLLIDING
        )

    def test_batch_entry_sanitizes_and_handles_duplicates(self):
        mapping = sanitize_names_uniquely(["Test Name (x)", "Test Name (x)", "a;b,c"])

        assert mapping["Test Name (x)"] == "Test_Name_x"
        assert mapping["a;b,c"] == "a_b_c"
        assert len(mapping) == 2

    def test_batch_entry_never_yields_empty_name(self):
        mapping = sanitize_names_uniquely(["   ", "(((("])
        assert all(name for name in mapping.values())
        assert len(set(mapping.values())) == 2


class TestMethodMismatch:
    """C-29：``add_entry`` 命中已有条目时不得静默丢弃本次传入的 method。"""

    @pytest.fixture
    def mapper(self):
        manager = NameMappingManager()
        manager.logger = RecordingLogger()
        return manager

    def _register_twice(self, mapper, on_method_mismatch=None):
        kwargs = (
            {}
            if on_method_mismatch is None
            else {"on_method_mismatch": on_method_mismatch}
        )
        first = mapper.add_entry("plain_name", NameShortenMode.ACCESSION)
        second = mapper.add_entry("plain_name", NameShortenMode.SERIAL, **kwargs)
        return first, second

    def test_warn_is_default_and_is_loud(self, mapper):
        first, second = self._register_twice(mapper)

        assert first == second == "S00001"
        assert mapper.logger.logged("WARNING", "add_entry 策略冲突")

    def test_raise_option(self, mapper):
        with pytest.raises(ValueError):
            self._register_twice(mapper, on_method_mismatch="raise")

    def test_regenerate_option_rebuilds_under_new_mode(self, mapper):
        mapper.add_entry("plain_name", NameShortenMode.SERIAL)
        rebuilt = mapper.add_entry(
            "plain_name", NameShortenMode.FULL, on_method_mismatch="regenerate"
        )

        assert rebuilt == "plain_name"
        # 旧短名让位：反查表只剩一条，不会留下指向已废弃条目的幽灵
        assert mapper.get_reverse_mapping_dict() == {"plain_name": "plain_name"}
        assert mapper.check_mapping_integrity() == []

    def test_same_mode_twice_stays_silent(self, mapper):
        mapper.add_entry("plain_name", NameShortenMode.FULL)
        mapper.add_entry("plain_name", NameShortenMode.FULL)

        assert not any(lvl == "WARNING" for lvl, _ in mapper.logger.records)

    def test_unknown_policy_rejected(self, mapper):
        mapper.add_entry("plain_name", NameShortenMode.FULL)
        with pytest.raises(ValueError):
            mapper.add_entry(
                "plain_name", NameShortenMode.SERIAL, on_method_mismatch="nonsense"
            )
