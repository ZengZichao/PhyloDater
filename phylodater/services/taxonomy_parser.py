"""
TaxonomyParser - GTDB/NCBI 分类学解析器

支持 GTDB 和 NCBI 两套分类学命名格式
支持从外部文件加载分类映射

格式A（嵌入式）：_d_Bacteria_p_Cyanobacteriota_...
格式B（表格分号式）：d__Archaea;p__Thermoproteota;...

加载契约（审阅报告 B-26）：``load_from_file`` **不会**在"整份表一行都没加载进来"
时返回一个看起来正常的 0。任何静默失效路径（表头缺 name 列、列名不被识别、
分类串全部解析失败、文件为空）都以 ``TaxonomyLoadError``（``ValueError`` 子类）
上抛；成功的加载同时返回条目数并在 ``last_load_report`` 留下逐行计数
（条目数 / 行数 / 各类跳过原因），调用方必须检查。
``check_label_consistency`` 在没有可用分类表时返回 ``None``——它绝不把全部叶
节点报成"已匹配"。
"""

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from ..infrastructure.logging import get_logger
from ..infrastructure.safe_io import safe_writer

#: 分类学解析结果：等级名 -> 该等级的取值。**缺失的等级以 None 显式入库**
#: （``phylodater check`` 的“格式B … species=None (缺失值处理正确)”就是在
#: 挨这条契约），所以值类型必须是 ``Optional[str]``；旧标注写成
#: ``Dict[str, str]`` 与实现不符，是十余条 assignment/var-annotated 报错的根源。
TaxonomyMap = Dict[str, Optional[str]]

# 默认分类级别映射（单字母前缀 -> 标准名称）
DEFAULT_TAXONOMY_LEVELS = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
    "k": "kingdom",
    "ss": "subspecies",
    "t": "strain",
}

# 反向映射：标准名称 -> 单字母前缀
DEFAULT_RANK_TO_PREFIX = {v: k for k, v in DEFAULT_TAXONOMY_LEVELS.items()}

# 多列表格里可充当"名称列"的列名别名（规范化后小写、去 BOM/前导 #）。
# 注意不能把 "species" 之类的等级列放进来——它是分类级别而不是行标识。
NAME_COLUMN_ALIASES: Tuple[str, ...] = (
    "name",
    "taxon",
    "taxon_name",
    "organism",
    "leaf",
    "leaf_name",
    "label",
    "seqname",
    "sequence_name",
    # GTDB metadata 表里真正的行标识（优先级低于上面这些，避免把
    # "taxon 写在 accession 之后"的表读错）
    "assembly_name",
    "accession",
)

# 两列表格里"第一列是表头而非数据"的识别词（第二列还解析不出分类串时才生效）
TWO_COLUMN_HEADER_TOKENS: Tuple[str, ...] = NAME_COLUMN_ALIASES + ("id", "#id")

# 存放整条分类串的列名（GTDB metadata 的 ``classification`` 等）
TAXONOMY_STRING_COLUMNS: Tuple[str, ...] = (
    "classification",
    "taxonomy",
    "lineage",
    "taxon_lineage",
    "full_taxonomy",
)


class TaxonomyLoadError(ValueError):
    """分类学表加载后没有任何可用条目（B-26：不得伪装成"加载成功，0 条"）。"""


@dataclass
class TaxonomyLoadReport:
    """``load_from_file`` 的结构化结果，供调用方做对账（B-26 / C-43）。

    ``entries`` 是**实际可用的去重条目数**（等于 ``len(_file_loaded_taxonomy)``），
    不再按行计数，因此重复行不会把它虚增。
    """

    file_path: Optional[Path] = None
    layout: str = "unknown"  # two_column | multi_column | unnamed_two_column
    header_detected: bool = False
    columns: List[str] = field(default_factory=list)
    name_column: Optional[str] = None
    rank_columns: List[str] = field(default_factory=list)
    string_column: Optional[str] = None
    ignored_columns: List[str] = field(default_factory=list)

    data_rows: int = 0  # 参与解析的非空行数（不含表头）
    entries: int = 0  # 去重后的可用条目数
    skipped_no_name: int = 0  # 名称列为空的行
    skipped_unparsed: int = 0  # 分类串解析不出任何级别
    skipped_empty_values: int = 0  # 只有空级别（如 's__'）的行
    skipped_malformed: int = 0  # 列数与表头不符的行
    short_rows: int = 0  # 少于表头列数的行（缺失单元格按空处理，不丢行）
    long_rows: int = 0  # 多于表头列数的行（多余单元格被忽略）
    duplicate_names: int = 0  # 与更早行重名的行（保留先入者）
    conflicting_duplicates: int = 0  # 其中取值还不一致的行
    failed_names: List[str] = field(default_factory=list)  # 供日志展示的样例

    @property
    def is_success(self) -> bool:
        """是否至少加载到一条可用条目。"""
        return self.entries > 0

    @property
    def skipped_total(self) -> int:
        return (
            self.skipped_no_name
            + self.skipped_unparsed
            + self.skipped_empty_values
            + self.skipped_malformed
            + self.duplicate_names
        )

    @property
    def has_disclosures(self) -> bool:
        """是否存在需要向用户披露的跳过/重复行。"""
        return self.skipped_total > 0

    def summary(self) -> str:
        """一行式对账文本：条目数 / 行数 / 各类跳过原因。"""
        parts = [
            f"layout={self.layout}",
            f"entries={self.entries}",
            f"data_rows={self.data_rows}",
        ]
        detail = []
        if self.skipped_unparsed:
            detail.append(f"unparsed={self.skipped_unparsed}")
        if self.skipped_empty_values:
            detail.append(f"empty_ranks_only={self.skipped_empty_values}")
        if self.skipped_no_name:
            detail.append(f"missing_name={self.skipped_no_name}")
        if self.skipped_malformed:
            detail.append(f"malformed_rows={self.skipped_malformed}")
        if self.duplicate_names:
            detail.append(
                f"duplicate_names={self.duplicate_names}"
                f"(conflicting={self.conflicting_duplicates})"
            )
        if detail:
            parts.append("skipped[" + ", ".join(detail) + "]")
        ragged = []
        if self.short_rows:
            ragged.append(f"short_rows={self.short_rows}")
        if self.long_rows:
            ragged.append(f"long_rows={self.long_rows}")
        if ragged:
            parts.append("ragged[" + ", ".join(ragged) + "]")
        if self.failed_names:
            sample = ", ".join(self.failed_names[:5])
            more = "…" if len(self.failed_names) > 5 else ""
            parts.append(f"first_failed=[{sample}{more}]")
        return " | ".join(parts)


@dataclass
class TaxonomyConsistencyReport:
    """``evaluate_label_consistency`` 的结果（B-26）。

    ``performed=False`` 表示这次检查**根本没有做**（没有可用分类表），
    此时三个清单都是空的，调用方必须把它读成"未核对"，而不是"全部匹配"。
    """

    performed: bool
    reason: str = ""
    tree_tips: int = 0
    table_entries: int = 0
    matched: List[str] = field(default_factory=list)
    unmatched_in_tree: List[str] = field(default_factory=list)
    unmatched_in_table: List[str] = field(default_factory=list)

    @property
    def fully_matched(self) -> bool:
        """检查已执行且两侧都没有未匹配项。"""
        return (
            self.performed
            and not self.unmatched_in_tree
            and not self.unmatched_in_table
        )

    def summary(self) -> str:
        if not self.performed:
            return f"一致性检查未执行（{self.reason}）"
        return (
            f"matched={len(self.matched)}/{self.tree_tips} tree tips, "
            f"{self.table_entries} table entries, "
            f"unmatched_in_tree={len(self.unmatched_in_tree)}, "
            f"unmatched_in_table={len(self.unmatched_in_table)}"
        )


class TaxonomyParser:
    """
    GTDB/NCBI 分类学解析器

    支持三种格式解析:
    1. 格式A（嵌入式）: GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_...
       分隔符格式: _{单字母前缀}_ （下划线数量固定）
    2. 格式B（表格分号式）: d__Bacteria;p__Cyanobacteriota;...
       分隔符格式: {单字母前缀}__
    3. NCBI 分号分隔格式: Bacteria;Proteobacteria;Gammaproteobacteria;...

    支持从外部 TSV 文件加载分类映射
    支持自定义分类级别前缀
    """

    def __init__(
        self,
        custom_levels: Optional[Dict[str, str]] = None,
        delimiter_mode: str = "reverse",
        source_priority: str = "table",
    ) -> None:
        """
        初始化分类学解析器

        Args:
            custom_levels: 自定义分类级别映射 {前缀: 级别名称}
                          例如: {'k': 'kingdom', 'ss': 'subspecies'}
            delimiter_mode: 解析模式 ('reverse', 'greedy', 'segment')
                - 'reverse': 从右向左解析，处理歧义前缀
                - 'greedy': 从左向右贪心解析
                - 'segment': 每个 _X_ 段严格对应一个级别
            source_priority: 分类信息来源优先级 ('table', 'embedded')
                - 'table': 文件分类优先
                - 'embedded': 名称中嵌入的分类优先
        """
        self._levels = {**DEFAULT_TAXONOMY_LEVELS}
        if custom_levels:
            self._levels.update(custom_levels)

        self._rank_to_prefix = {v: k for k, v in self._levels.items()}

        self._prefix_chars = "".join(sorted(set(self._levels.keys())))

        self._compile_patterns()

        self._cache: Dict[str, TaxonomyMap] = {}
        self._file_loaded_taxonomy: Dict[str, TaxonomyMap] = {}
        self._unmatched_count = 0
        self._warning_threshold = 5
        self._warning_issued = False
        self._loaded_file_path: Optional[Path] = None
        self._taxonomy_source: Dict[str, str] = {}
        self._delimiter_mode = delimiter_mode
        self._source_priority = source_priority
        # 最近一次 load_from_file 的结构化结果（B-26）；None 表示从未尝试加载
        self.last_load_report: Optional[TaxonomyLoadReport] = None
        self.logger = get_logger()

    def _compile_patterns(self) -> None:
        """编译正则表达式模式"""
        self.UNDERSCORE_PATTERN = re.compile(r"_([a-z]+)_" r"((?:(?!_[a-z]+_).)+)")

        # segment 模式专用分隔符：前缀必须是「已知分类级别前缀」（来自 self._levels）。
        # 这样 "_d_X_Y_" 中 Y 不是级别前缀时不会被误判为新的级别分隔符，
        # 从而保证含下划线的级值（如 "X_Y"）不被静默丢级。
        # 按前缀长度降序排列，确保 "ss" 优先于 "s" 等更长前缀被先匹配。
        seg_prefixes = "|".join(
            re.escape(p) for p in sorted(self._levels, key=len, reverse=True)
        )
        self.UNDERSCORE_SEGMENT_PATTERN = re.compile(
            r"_(" + seg_prefixes + r")_((?:(?!_(" + seg_prefixes + r")_).)+)"
        )

        self.DOUBLE_UNDERSCORE_PATTERN = re.compile(r"([a-z]{1,2})__([^;]*)")

        self.NCBI_SEMICOLON_PATTERN = re.compile(
            r"^(?P<domain>[^;]+)"
            r"(?:;(?P<phylum>[^;]+))?"
            r"(?:;(?P<class>[^;]+))?"
            r"(?:;(?P<order>[^;]+))?"
            r"(?:;(?P<family>[^;]+))?"
            r"(?:;(?P<genus>[^;]+))?"
            r"(?:;(?P<species>[^;]+))?"
            r"(?:;(?P<strain>[^;]+))?$"
        )

        self.NCBI_TAXID_PATTERN = re.compile(r"taxid:?(\d+)", re.IGNORECASE)

        self.NCBI_ACCESSION_PATTERN = re.compile(
            r"([A-Z]{3}[-_]?\d{5,}(?:\.\d+)?)", re.IGNORECASE
        )

    def clear_cache(self) -> None:
        """清除解析缓存"""
        self._cache.clear()
        self._unmatched_count = 0
        self._warning_issued = False

    def _store_rank(
        self,
        result: TaxonomyMap,
        name: str,
        prefix: str,
        value: str,
        conflicts: Optional[Dict[str, List[Optional[str]]]] = None,
    ) -> bool:
        """将 prefix→rank_name 解析结果写入 result（统一空值处理）。

        供各 ``_parse_underscore_*`` 与 ``_parse_double_underscore_format`` 复用，
        消除重复的 ``if prefix in self._levels: rank_name = ...`` 赋值片段
        （#16 可维护性收敛）。返回 True 表示 prefix 属于已知级别并已写入。

        同一等级前缀在名称里出现两次时（如 ``…_g_GenusA_g_GenusB…``）后值覆盖前值，
        但**不再静默**覆盖：把该级别的取值序列记进 ``conflicts``，由调用方在解析
        出口统一发 WARNING（C-33）。``conflicts`` 为 None 时只覆盖不留痕。
        """
        if prefix not in self._levels:
            return False
        rank_name = self._levels[prefix]
        new_value = value if value else None
        if rank_name in result and result[rank_name] != new_value:
            if conflicts is not None:
                seen = conflicts.setdefault(rank_name, [result[rank_name]])
                seen.append(new_value)
        if new_value is not None:
            result[rank_name] = new_value
        else:
            result[rank_name] = None
            self.logger.debug(f"Empty value for rank '{rank_name}' in: {name}")
        return True

    def _warn_rank_conflicts(
        self, source: str, conflicts: Dict[str, List[Optional[str]]]
    ) -> None:
        """把 ``_store_rank`` 收集的重复等级前缀冲突一次性披露出去（C-33）。"""
        for rank_name, values in conflicts.items():
            self.logger.warning(
                f"分类串中同一等级前缀重复出现 for '{source}': rank '{rank_name}' "
                f"依次取值为 {values}，最终保留最后出现的 "
                f"{values[-1]!r}。请核对名称/表格是否混入了两个同名级（"
                f"如 '_g_GenusA_g_GenusB'），必要时改用 --taxonomy-delimiter-mode "
                f"segment 或提供外部分类表。"
            )

    def parse(self, name: str) -> Optional[TaxonomyMap]:
        """
        解析叶节点名称中的分类学信息

        优先级（受 source_priority 控制）:
        - source_priority='table': 文件分类 > 缓存 > 格式B > 格式A > NCBI
        - source_priority='embedded': 缓存 > 格式B > 格式A > NCBI > 文件分类

        当同一节点同时具有嵌入式和表格两种来源时，按优先级选择主结果，
        但如有冲突（重叠级别值不同）发出 WARNING。

        Args:
            name: 叶节点名称

        Returns:
            分类学字典，解析失败返回 None
        """
        if not name:
            return None

        if name in self._cache:
            return self._cache[name]

        has_table = name in self._file_loaded_taxonomy

        if self._source_priority == "table":
            if has_table:
                result = self._file_loaded_taxonomy[name]
                # 冲突检测：同时存在表格和嵌入式来源
                try:
                    embedded_result = self._parse_embedded(name)
                    if embedded_result:
                        self._check_taxonomy_conflict(
                            name, embedded_result, result, "table"
                        )
                except ValueError:
                    pass  # 嵌入式解析失败，不影响表格结果
                self._cache[name] = result
                self._taxonomy_source[name] = "file"
                return result

            embedded_result = self._parse_embedded(name)
            if embedded_result:
                self._cache[name] = embedded_result
                return embedded_result
        else:
            embedded_result = self._parse_embedded(name)
            if embedded_result:
                # 冲突检测：同时存在嵌入式和表格来源
                if has_table:
                    table_result = self._file_loaded_taxonomy[name]
                    self._check_taxonomy_conflict(
                        name, embedded_result, table_result, "embedded"
                    )
                self._cache[name] = embedded_result
                return embedded_result

            if has_table:
                result = self._file_loaded_taxonomy[name]
                self._cache[name] = result
                self._taxonomy_source[name] = "file"
                return result

        self._unmatched_count += 1
        if (
            self._unmatched_count >= self._warning_threshold
            and not self._warning_issued
        ):
            self.logger.warning(
                f"{self._unmatched_count} taxa could not be parsed. "
                "Consider providing external taxonomy file via --taxonomy-file."
            )
            self._warning_issued = True

        return None

    def _check_taxonomy_conflict(
        self, name: str, embedded: TaxonomyMap, table: TaxonomyMap, priority: str
    ) -> None:
        """
        检测嵌入式和表格来源的分类学信息冲突

        比对两来源中重叠的级别，如有任何级别的值不同则发出 WARNING。
        不影响最终返回的结果（仍按优先级选择）。

        Args:
            name: 叶节点名称
            embedded: 嵌入式解析结果
            table: 表格解析结果
            priority: 主来源 ('table' 或 'embedded')
        """
        common_levels = set(embedded.keys()) & set(table.keys())
        for level in common_levels:
            embedded_val = embedded.get(level)
            table_val = table.get(level)
            if embedded_val != table_val:
                self.logger.warning(
                    f"分类学来源冲突 for '{name}': level '{level}' "
                    f"embedded='{embedded_val}' vs table='{table_val}'，"
                    f"使用 {priority} 来源"
                )

    def _parse_embedded(self, name: str) -> Optional[TaxonomyMap]:
        # 重复等级前缀的冲突在这里统一收集并在出口披露（C-33），
        # 使嵌入式解析的三条分支共用同一套留痕规则。
        conflicts: Dict[str, List[Optional[str]]] = {}

        result = self._parse_double_underscore_format(name, conflicts=conflicts)
        if result:
            self._taxonomy_source[name] = "gtdb_double_underscore"
            self._warn_rank_conflicts(name, conflicts)
            return result

        result = self._parse_underscore_format(name)
        if result:
            self._taxonomy_source[name] = "gtdb_underscore"
            return result

        result = self._parse_ncbi_semicolon_format(name)
        if result:
            self._taxonomy_source[name] = "ncbi"
            return result

        self._warn_rank_conflicts(name, conflicts)
        return None

    def _parse_underscore_format(self, name: str) -> Optional[TaxonomyMap]:
        """
        解析格式A（嵌入式）: _d_Bacteria_p_Cyanobacteriota_...

        支持三种模式：
        - 'reverse': 从右向左解析，处理歧义前缀
        - 'greedy': 从左向右贪心解析
        - 'segment': 每个 _X_ 段严格对应一个级别
        """
        conflicts: Dict[str, List[Optional[str]]] = {}
        if self._delimiter_mode == "segment":
            result = self._parse_underscore_segment(name, conflicts)
        elif self._delimiter_mode == "reverse":
            result = self._parse_underscore_reverse(name, conflicts)
        else:
            result = self._parse_underscore_greedy(name, conflicts)

        if result and conflicts:
            # C-33：重复等级前缀的"后值覆盖前值"必须有痕迹
            self._warn_rank_conflicts(name, conflicts)
        return result

    def _parse_underscore_greedy(
        self, name: str, conflicts: Optional[Dict[str, List[Optional[str]]]] = None
    ) -> Optional[TaxonomyMap]:
        result: TaxonomyMap = {}

        for match in self.UNDERSCORE_PATTERN.finditer(name):
            prefix = match.group(1)
            value = match.group(2).strip("_").strip()

            self._store_rank(result, name, prefix, value, conflicts)

        return result if result else None

    def _parse_underscore_reverse(
        self, name: str, conflicts: Optional[Dict[str, List[Optional[str]]]] = None
    ) -> Optional[TaxonomyMap]:
        greedy_result = self._parse_underscore_greedy(name, conflicts)
        if greedy_result is None:
            return None

        matches = list(self.UNDERSCORE_PATTERN.finditer(name))
        ambiguous = False
        for match in matches:
            prefix = match.group(1)
            if prefix in self._levels:
                for other_prefix, other_rank in self._levels.items():
                    if other_prefix != prefix and other_prefix.startswith(prefix):
                        ambiguous = True
                        break
            if ambiguous:
                break

        if not ambiguous:
            return greedy_result

        reverse_result: TaxonomyMap = {}
        # 与 greedy 走同一个写入出口，避免两条路径的覆盖/空值语义漂移（C-33）
        reverse_conflicts: Dict[str, List[Optional[str]]] = {}
        for match in reversed(matches):
            prefix = match.group(1)
            value = match.group(2).strip("_").strip()

            self._store_rank(reverse_result, name, prefix, value, reverse_conflicts)

        if reverse_result != greedy_result:
            self.logger.warning(
                f"Ambiguous parsing detected for '{name}': "
                f"greedy={greedy_result}, reverse={reverse_result}. "
                f"Using reverse parse result."
            )
            return reverse_result

        self.logger.warning(
            f"Ambiguous prefix detected in '{name}', "
            f"but both parse directions agree."
        )
        return greedy_result

    def _parse_underscore_segment(
        self, name: str, conflicts: Optional[Dict[str, List[Optional[str]]]] = None
    ) -> Optional[TaxonomyMap]:
        result: TaxonomyMap = {}

        for match in self.UNDERSCORE_SEGMENT_PATTERN.finditer(name):
            prefix = match.group(1)
            # 与 greedy 模式保持一致：去掉值两端可能残留的分隔符下划线，
            # 例如 "_d_X_Y_" 的级值应解析为 "X_Y" 而非 "X_Y_"。
            value = match.group(2).strip("_").strip()

            self._store_rank(result, name, prefix, value, conflicts)

        return result if result else None

    def _parse_double_underscore_format(
        self,
        name: str,
        taxonomy_sep: str = ";",
        conflicts: Optional[Dict[str, List[Optional[str]]]] = None,
    ) -> Optional[TaxonomyMap]:
        """
        解析格式B（表格分号式）: d__Bacteria;p__Thermoproteota;...

        特征：
        - 分隔符格式: {单字母前缀}__
        - 各级别用 taxonomy_sep 分隔（默认分号）
        - 种水平可缺失值（如 s__ 解析为 None）
        """
        if "__" not in name:
            return None

        result: TaxonomyMap = {}

        if taxonomy_sep == ";":
            pattern = self.DOUBLE_UNDERSCORE_PATTERN
        else:
            sep_pattern = re.escape(taxonomy_sep)
            pattern = re.compile(r"([a-z]{1,2})__([^" + sep_pattern + r"]*)")

        for match in pattern.finditer(name):
            prefix = match.group(1)
            value = match.group(2).strip()

            for i, ch in enumerate(value):
                if ord(ch) < 0x20 and ch not in ("\t", "\n"):
                    raise ValueError(
                        f"Illegal control character (0x{ord(ch):02x}) in taxonomy value "
                        f"for prefix '{prefix}' in: {name}"
                    )

            if prefix in self._levels:
                # 净化检查：值不得包含 __ 或分隔符
                if "__" in value:
                    self.logger.error(
                        f"Malformed taxonomy value contains '__' for prefix '{prefix}' in: {name}"
                    )
                    continue
                if taxonomy_sep in value:
                    self.logger.error(
                        f"Malformed taxonomy value contains separator '{taxonomy_sep}' "
                        f"for prefix '{prefix}' in: {name}"
                    )
                    continue

                self._store_rank(result, name, prefix, value, conflicts)

        return result if result else None

    def _parse_ncbi_semicolon_format(self, name: str) -> Optional[TaxonomyMap]:
        """
        解析 NCBI 分号分隔格式: Bacteria;Proteobacteria;...

        特征：
        - 无前缀标记
        - 各级别用分号分隔
        - 按顺序映射到 domain, phylum, class, ...
        """
        if ";" not in name:
            return None

        parts = name.split(";")
        if len(parts) < 2:
            return None

        ranks = [
            "domain",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "strain",
        ]
        result: TaxonomyMap = {}

        for i, part in enumerate(parts):
            part = part.strip()
            if part and i < len(ranks):
                result[ranks[i]] = part

        domain = result.get("domain")
        if domain in ("Bacteria", "Archaea", "Eukarya", "Eukaryota"):
            normalized = domain.capitalize()
            if normalized == "Eukaryota":
                normalized = "Eukarya"
            result["domain"] = normalized

        return result if result else None

    def get_taxa_by_domain(self, names: List[str], domain: str) -> List[str]:
        """
        获取指定域的所有叶节点

        Args:
            names: 叶节点名称列表
            domain: 域名称（Bacteria, Archaea）

        Returns:
            匹配的叶节点名称列表
        """
        result = []
        domain_lower = domain.lower()

        for name in names:
            taxonomy = self.parse(name)
            # ``domain`` 这个等级可以显式存 None（缺失级），所以取完再归一成空串。
            if taxonomy and (taxonomy.get("domain") or "").lower() == domain_lower:
                result.append(name)

        return result

    def get_taxa_by_classification(
        self, names: List[str], rank: str, value: str
    ) -> List[str]:
        """
        根据分类级别获取叶节点

        Args:
            names: 叶节点名称列表
            rank: 分类级别
            value: 分类值

        Returns:
            匹配的叶节点名称列表
        """
        result = []
        value_lower = value.lower()

        for name in names:
            taxonomy = self.parse(name)
            if taxonomy:
                rank_value = taxonomy.get(rank)
                if rank_value is not None and rank_value.lower() == value_lower:
                    result.append(name)

        return result

    def get_taxa_with_missing_rank(self, names: List[str], rank: str) -> List[str]:
        """获取指定分类级别缺失的叶节点"""
        result = []
        for name in names:
            taxonomy = self.parse(name)
            if taxonomy and (rank not in taxonomy or taxonomy[rank] is None):
                result.append(name)
        return result

    def get_taxonomy_statistics(self, names: List[str]) -> Dict[str, Dict[str, int]]:
        """获取分类学统计信息"""
        stats = {}
        for name in names:
            taxonomy = self.parse(name)
            if not taxonomy:
                continue
            for rank, value in taxonomy.items():
                if rank not in stats:
                    stats[rank] = {"total": 0, "filled": 0, "missing": 0}
                stats[rank]["total"] += 1
                if value is not None and value:
                    stats[rank]["filled"] += 1
                else:
                    stats[rank]["missing"] += 1
        return stats

    def is_target_in_name(self, name: str, target: str) -> bool:
        """检查目标是否在名称中（词边界匹配）"""
        pattern = rf"(_|^){re.escape(target)}(_|$)"
        return bool(re.search(pattern, name, re.IGNORECASE))

    def is_exact_match(self, name: str, target: str) -> bool:
        """精确分词匹配"""
        normalized = re.sub(r"[;_.\-\s]+", " ", name.lower())
        name_tokens = set(normalized.split())
        target_normalized = target.lower()
        if " " not in target_normalized:
            return target_normalized in name_tokens
        pattern = rf"\b{re.escape(target_normalized)}\b"
        return bool(re.search(pattern, normalized))

    def evaluate_label_consistency(
        self, tree_tip_names: Sequence[str]
    ) -> TaxonomyConsistencyReport:
        """分类表 ↔ 树标签一致性检查的结构化入口（B-26）。

        与 :meth:`check_label_consistency` 的唯一区别是返回值永远不会把
        "没做检查"伪装成"全部匹配"：检查未执行时 ``performed=False`` 且
        ``reason`` 说明原因。

        Args:
            tree_tip_names: 树的叶节点名称列表

        Returns:
            ``TaxonomyConsistencyReport``
        """
        tips = list(tree_tip_names)

        if not self._file_loaded_taxonomy:
            reason = (
                "没有可用的外部分类表"
                if self._loaded_file_path is None
                else f"分类表 {self._loaded_file_path} 未提供任何条目"
            )
            self.logger.warning(
                f"Taxonomy consistency check NOT PERFORMED: {reason}. "
                f"{len(tips)} tree tip(s) remain UNVERIFIED against any taxonomy "
                f"table —— 这不是'全部匹配'。"
            )
            return TaxonomyConsistencyReport(
                performed=False,
                reason=reason,
                tree_tips=len(tips),
                table_entries=0,
            )

        tree_set = set(tips)
        table_set = set(self._file_loaded_taxonomy.keys())

        matched = sorted(tree_set & table_set)
        unmatched_in_tree = sorted(tree_set - table_set)
        unmatched_in_table = sorted(table_set - tree_set)

        if unmatched_in_tree:
            self.logger.warning(
                f"{len(unmatched_in_tree)} tree tip(s) not found in taxonomy table: "
                f"{unmatched_in_tree[:10]}{'...' if len(unmatched_in_tree) > 10 else ''}"
            )
        if unmatched_in_table:
            self.logger.warning(
                f"{len(unmatched_in_table)} taxonomy table entry(ies) not found in tree: "
                f"{unmatched_in_table[:10]}{'...' if len(unmatched_in_table) > 10 else ''}"
            )
        if not matched:
            self.logger.warning(
                f"0 of {len(tree_set)} tree tip(s) match the {len(table_set)} "
                f"taxonomy table entry(ies) —— 分类表与树可能来自不同数据集。"
            )

        return TaxonomyConsistencyReport(
            performed=True,
            tree_tips=len(tree_set),
            table_entries=len(table_set),
            matched=matched,
            unmatched_in_tree=unmatched_in_tree,
            unmatched_in_table=unmatched_in_table,
        )

    def check_label_consistency(
        self, tree_tip_names: List[str]
    ) -> Optional[Tuple[List[str], List[str], List[str]]]:
        """
        检查分类表与树标签的一致性

        Args:
            tree_tip_names: 树的叶节点名称列表

        Returns:
            ``(matched_names, unmatched_in_tree, unmatched_in_table)``，三个清单
            均已排序（C-45：``list(set)`` 的顺序随 PYTHONHASHSEED 变化，不可复现）；
            **没有可用分类表时返回 ``None``**，表示检查未执行（B-26）。

        Note:
            旧实现在分类表为空时返回 ``(全部叶名, [], [])``，于是"表↔树一致性"
            这项质量检查在表完全失效时呈现为 100% 通过。调用方现在必须显式处理
            ``None``（打印"未做一致性检查"并以 ``EXIT_DATA_ERROR`` 之类退出码
            终止或降级），不得把它当成"全部匹配"。
        """
        report = self.evaluate_label_consistency(tree_tip_names)
        if not report.performed:
            return None
        return report.matched, report.unmatched_in_tree, report.unmatched_in_table

    def load_from_file(
        self,
        file_path: Path,
        table_sep: Optional[str] = None,
        taxonomy_sep: str = ";",
    ) -> int:
        """
        从文件加载分类映射

        支持三种格式：
        1. 两列格式（制表符或逗号分隔）：
           name<TAB>d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...
           或
           name,d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...

        2. 多列 TSV 格式：
           name<TAB>domain<TAB>phylum<TAB>class<TAB>order<TAB>family<TAB>genus<TAB>species
           列名不区分大小写、允许前后空格/引号/前导 ``#``；名称列也可写作
           ``taxon`` / ``organism`` / ``label`` 等（见 ``NAME_COLUMN_ALIASES``）；
           ``kingdom`` / ``subspecies`` / ``strain`` 以及整条分类串列
           ``classification`` / ``taxonomy``（GTDB metadata 写法）同样会被读取。

        3. 双下划线格式的分类字符串（如 GTDB 格式）：
           d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__

        Args:
            file_path: 分类映射文件路径
            table_sep: 表格列分隔符，None 时自动检测
            taxonomy_sep: 分类学字符串内部的级别分隔符（默认 ;，可设为 | 等）

        Returns:
            加载的分类条目数量（**去重后的条目数**，不是行数；同一数字也保存在
            ``self.last_load_report.entries``，报告里还有逐类跳过原因）

        Raises:
            FileNotFoundError: 文件不存在
            TaxonomyLoadError: 加载后没有任何可用条目（``ValueError`` 的子类，
                因此只捕获 ``ValueError`` 的既有调用方同样能接住）。覆盖：文件为空、
                表头没有名称列、没有任何可识别的分类列、所有行的分类串都解析失败。
                旧实现对这些情形只打一条 ``Loaded 0 taxonomy entries`` 的 INFO 就
                返回 0（B-26 的实测成因）。
            ValueError: 文件格式错误 / 含非法控制字符

        Note:
            全部校验通过后才提交内部状态；失败时既有分类表保持不变，避免
            "加载坏文件顺手清空好表"的连带失效。提交时同时清空解析缓存，
            因为表格来源的优先级高于名称内嵌来源。
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"Taxonomy file not found: {file_path}")

        report = TaxonomyLoadReport(file_path=file_path)
        new_table: Dict[str, TaxonomyMap] = {}

        try:
            with open(file_path, "r", newline="", encoding="utf-8-sig") as f:
                lines = [line.strip() for line in f if line.strip()]

            if not lines:
                raise TaxonomyLoadError(f"Taxonomy file is empty: {file_path}")

            delimiter = (
                table_sep
                if table_sep is not None
                else ("\t" if "\t" in lines[0] else ",")
            )

            first_line_parts = lines[0].split(delimiter)
            columns = self._normalize_columns(first_line_parts)
            report.columns = columns

            known_ranks = self._rank_column_names()
            rank_columns = [col for col in columns if col in known_ranks]
            string_column = next(
                (col for col in columns if col in TAXONOMY_STRING_COLUMNS), None
            )
            name_column = self._find_name_column(columns)
            is_two_column = len(first_line_parts) == 2

            if is_two_column and not rank_columns:
                # 两列 = 名称 + 整条分类串（第二列名可以是 taxonomy / classification，
                # 由表头词识别后被当作名称列的别名，不参与两列判定）
                report.layout = "two_column"
                self.logger.info("Detected two-column format with taxonomy strings")
                self._load_two_column_rows(
                    lines, delimiter, taxonomy_sep, new_table, report
                )
            elif name_column or rank_columns or string_column:
                # 表头被识别出来：必需列缺失即失败，绝不逐行静默 skip（B-26）
                report.layout = "multi_column"
                report.header_detected = True
                if name_column is None:
                    raise TaxonomyLoadError(
                        f"分类表 {file_path} 的表头没有可用的名称列，无法把条目对应到树标签。"
                        f" 实际列: {columns}；已识别到的分类级别列: {rank_columns or '无'}。"
                        f" 名称列需命名为 {' / '.join(NAME_COLUMN_ALIASES[:5])} 之一"
                        f"（大小写、前后空格与前导 # 不限）。"
                    )
                if not rank_columns and string_column is None:
                    raise TaxonomyLoadError(
                        f"分类表 {file_path} 有名称列 {name_column!r}，但没有任何可识别的"
                        f"分类级别列。 实际列: {columns}；可用列名: "
                        f"{sorted(known_ranks)}，或整条分类串列 "
                        f"{list(TAXONOMY_STRING_COLUMNS)}；自定义级别请用 "
                        f"custom_levels / --taxonomy-levels。"
                    )
                report.name_column = name_column
                report.rank_columns = list(rank_columns)
                report.string_column = string_column
                report.ignored_columns = [
                    col
                    for col in columns
                    if col not in known_ranks
                    and col != name_column
                    and col != string_column
                ]
                if report.ignored_columns:
                    self.logger.debug(
                        "Taxonomy table columns ignored (not taxonomy ranks): "
                        f"{report.ignored_columns}"
                    )
                self._load_multi_column_rows(
                    lines,
                    delimiter,
                    taxonomy_sep,
                    new_table,
                    report,
                    name_column,
                    rank_columns,
                    string_column,
                )
            else:
                # 无表头：退回"第一列名字 + 其余部分当作分类串"的宽松解析
                report.layout = "unnamed_two_column"
                self.logger.info(
                    "No taxonomy header detected; parsing as name + taxonomy string rows"
                )
                self._load_two_column_rows(
                    lines, delimiter, taxonomy_sep, new_table, report
                )

            report.entries = len(new_table)
        except TaxonomyLoadError:
            self.last_load_report = report
            raise
        except csv.Error as e:
            self.last_load_report = report
            raise ValueError(f"Invalid taxonomy file format: {e}") from e
        except Exception as e:
            self.last_load_report = report
            raise ValueError(f"Failed to parse taxonomy file: {e}") from e

        self.last_load_report = report

        if report.entries == 0:
            # B-26：空/失败的加载绝不再以"成功返回 0"的形式交给调用方
            raise TaxonomyLoadError(
                f"分类学表未加载到任何可用条目: {file_path}\n"
                f"  {report.summary()}\n"
                f"  实际列: {report.columns}\n"
                "  可能成因: (1) 表头缺少名称列 (name/taxon/…); "
                "(2) 等级列名不被识别; "
                "(3) 两列格式的第二列不是 d__X;p__Y 形式的分类串; "
                "(4) 文件只有表头、没有数据行。\n"
                "  调用方必须终止本次运行或显式降级：分类表缺失时 "
                "check_label_consistency() 返回 None（表示检查未执行），"
                "绝不能把全部叶节点当成已核对。"
            )

        # 提交（失败路径不会走到这里，既有分类表因此得以保留）
        self._file_loaded_taxonomy = new_table
        self._loaded_file_path = file_path
        self.clear_cache()  # 表已换，旧的解析结果不能继续命中缓存

        self.logger.info(
            f"Loaded {report.entries} taxonomy entries from {report.data_rows} "
            f"data row(s) in {file_path} [{report.summary()}]"
        )
        if report.has_disclosures:
            self.logger.warning(
                f"{report.skipped_total} of {report.data_rows} taxonomy table row(s) "
                f"were NOT used [{report.summary()}]. Check the header and the taxonomy "
                f"strings: silently ignoring rows would make rank-based filtering and "
                f"the table/tree consistency check under-report coverage."
            )
        return report.entries

    @staticmethod
    def _normalize_column(token: str) -> str:
        """规范化单个表头 token：去 BOM、前后空白、引号与前导 ``#``，再小写。"""
        cleaned = token.replace("\ufeff", "").strip().strip("\"'").lstrip("#").strip()
        return cleaned.lower()

    @classmethod
    def _normalize_columns(cls, tokens: Sequence[str]) -> List[str]:
        """规范化一整行列名。

        旧实现判定表头用"整行小写后取子串"、取值却用精确键名，两个口径不一致，
        于是 ``Name`` / `` Domain`` / ``#taxon``（GTDB metadata 的真实写法）这类
        表头会被识别成"有表头"却一列都取不到值——正是 B-26 的实测成因。
        """
        return [cls._normalize_column(token) for token in tokens]

    def _find_name_column(self, columns: Sequence[str]) -> Optional[str]:
        """在规范化后的列名里挑出名称列（按 ``NAME_COLUMN_ALIASES`` 的优先序）。"""
        present = set(columns)
        for alias in NAME_COLUMN_ALIASES:
            if alias in present:
                return alias
        return None

    def _rank_column_names(self) -> Set[str]:
        """多列表格里可识别的分类级别列名（来自当前生效的级别映射）。"""
        return set(self._levels.values())

    def _record_entry(
        self,
        table: Dict[str, TaxonomyMap],
        name: str,
        taxonomy: Optional[TaxonomyMap],
        report: TaxonomyLoadReport,
        where: str,
    ) -> bool:
        """把一行结果并入待提交的分类表，逐行留痕（B-26 / C-43）。

        三种"看起来加载了其实没用"的情形分别计数：

        - 名称列为空 → ``skipped_no_name``
        - 分类串解析不出任何级别 → ``skipped_unparsed``
        - 只有空级别（GTDB 常见的 ``s__`` 全空行）→ ``skipped_empty_values``；
          旧实现把这种行也算进"已加载条目"，于是 ``Loaded N`` 既高估又低估
        - 与更早的行重名 → 保留先入者并计 ``duplicate_names``，取值还不一致时
          再计 ``conflicting_duplicates`` 并 WARNING（C-43 的"后写者静默胜出"）
        """
        if not name:
            report.skipped_no_name += 1
            return False
        if not taxonomy:
            report.skipped_unparsed += 1
            if name not in report.failed_names:
                report.failed_names.append(name)
            return False
        if not any(value not in (None, "") for value in taxonomy.values()):
            report.skipped_empty_values += 1
            if name not in report.failed_names:
                report.failed_names.append(name)
            return False

        previous = table.get(name)
        if previous is not None:
            report.duplicate_names += 1
            if previous != taxonomy:
                report.conflicting_duplicates += 1
                self.logger.warning(
                    f"Duplicate taxonomy entry for '{name}' at {where}: kept the "
                    f"earlier {previous}, ignored {taxonomy} "
                    f"(C-43: duplicate rows no longer overwrite silently)"
                )
            else:
                self.logger.debug(
                    f"Identical duplicate taxonomy row for '{name}' at {where}: merged"
                )
            return False

        table[name] = taxonomy
        return True

    def _load_two_column_rows(
        self,
        lines: List[str],
        delimiter: str,
        taxonomy_sep: str,
        table: Dict[str, TaxonomyMap],
        report: TaxonomyLoadReport,
    ) -> None:
        """两列格式（name + 分类串）逐行加载。"""
        for line_no, line in enumerate(lines, start=1):
            parts = line.split(delimiter, 1)
            if len(parts) != 2:
                report.data_rows += 1
                report.skipped_malformed += 1
                self.logger.warning(f"Skipping malformed line: {line[:50]}...")
                continue

            name = parts[0].strip()
            taxonomy_str = parts[1].strip()
            taxonomy = (
                self._parse_taxonomy_string(taxonomy_str, taxonomy_sep=taxonomy_sep)
                if taxonomy_str
                else None
            )

            if (
                line_no == 1
                and not taxonomy
                and self._normalize_column(name) in TWO_COLUMN_HEADER_TOKENS
            ):
                # 表头行（如 "name<TAB>taxonomy"）：既不计入数据行，也不报成解析失败
                report.header_detected = True
                report.name_column = self._normalize_column(name)
                self.logger.debug(f"Two-column header detected: {line!r}")
                continue

            report.data_rows += 1
            self._record_entry(table, name, taxonomy, report, f"line {line_no}")

    def _load_multi_column_rows(
        self,
        lines: List[str],
        delimiter: str,
        taxonomy_sep: str,
        table: Dict[str, TaxonomyMap],
        report: TaxonomyLoadReport,
        name_column: str,
        rank_columns: List[str],
        string_column: Optional[str],
    ) -> None:
        """多列带表头格式加载（按列名取值，不依赖列序）。"""
        rows = list(csv.reader(lines, delimiter=delimiter))
        header = [self._normalize_column(col) for col in rows[0]]
        index_of: Dict[str, int] = {}
        for i, col in enumerate(header):
            index_of.setdefault(col, i)

        name_idx = index_of[name_column]
        rank_idx = [(index_of[rank], rank) for rank in rank_columns if rank in index_of]
        string_idx = index_of.get(string_column) if string_column else None

        for line_no, fields in enumerate(rows[1:], start=2):
            report.data_rows += 1

            if len(fields) < len(header):
                # 少于表头列数：缺失单元格按空处理（与旧实现的 csv.DictReader 同义，
                # 尾部空列在真实 TSV 里极其常见，直接丢行会让分类单元静默消失）。
                report.short_rows += 1
                fields = fields + [""] * (len(header) - len(fields))
            elif len(fields) > len(header):
                report.long_rows += 1
                fields = fields[: len(header)]

            name = fields[name_idx].strip()

            taxonomy: TaxonomyMap = {}
            if string_idx is not None:
                parsed = self._parse_taxonomy_string(
                    fields[string_idx].strip(), taxonomy_sep=taxonomy_sep
                )
                if parsed:
                    taxonomy.update(parsed)
            for idx, rank in rank_idx:
                value = fields[idx].strip()
                if value:
                    taxonomy[rank] = value

            self._record_entry(table, name, taxonomy, report, f"line {line_no}")

        if report.short_rows or report.long_rows:
            self.logger.warning(
                f"{report.short_rows} row(s) had fewer columns than the header "
                f"(missing cells read as empty) and {report.long_rows} row(s) had "
                f"more (extra cells ignored). Header: {header}"
            )

    def _parse_taxonomy_string(
        self, taxonomy_str: str, taxonomy_sep: str = ";"
    ) -> Optional[TaxonomyMap]:
        """
        解析分类学字符串

        支持格式：
        - 格式B: d__Bacteria;p__Proteobacteria;...（双下划线）
        - 格式A: _d_Bacteria_p_Proteobacteria_...（嵌入式，带分号分隔）
        - 格式B变体: d_Bacteria;p_Proteobacteria;...（单下划线）

        Args:
            taxonomy_str: 分类学字符串
            taxonomy_sep: 分类学字符串内部的级别分隔符（默认 ;）

        Returns:
            分类学字典，解析失败返回 None
        """
        if not taxonomy_str:
            return None

        conflicts: Dict[str, List[Optional[str]]] = {}
        result = self._parse_double_underscore_format(
            taxonomy_str, taxonomy_sep=taxonomy_sep, conflicts=conflicts
        )
        if result:
            if conflicts:
                self._warn_rank_conflicts(taxonomy_str, conflicts)
            return result

        if taxonomy_sep in taxonomy_str:
            parts = taxonomy_str.split(taxonomy_sep)
            taxonomy: TaxonomyMap = {}

            for part in parts:
                part = part.strip()
                if not part:
                    continue

                for prefix, rank_name in self._levels.items():
                    double_pattern = f"{prefix}__"
                    single_pattern = f"{prefix}_"

                    if part.startswith(double_pattern):
                        value = part[len(double_pattern) :].strip()
                        self._note_level_value(
                            taxonomy, rank_name, value if value else None, conflicts
                        )
                        break
                    elif part.startswith(single_pattern):
                        value = part[len(single_pattern) :].strip()
                        self._note_level_value(
                            taxonomy, rank_name, value if value else None, conflicts
                        )
                        break

            if conflicts:
                self._warn_rank_conflicts(taxonomy_str, conflicts)
            return taxonomy if taxonomy else None

        return None

    def _note_level_value(
        self,
        taxonomy: TaxonomyMap,
        rank_name: str,
        value: Optional[str],
        conflicts: Dict[str, List[Optional[str]]],
    ) -> None:
        """写入一个级别值，重复级别记入 ``conflicts``（C-33 的表格串侧）。

        与 :meth:`_store_rank` 保持同一套"后值覆盖前值 + 留痕"语义，避免
        ``_parse_taxonomy_string`` 与嵌入式解析两条路径行为漂移。
        """
        if rank_name in taxonomy and taxonomy[rank_name] != value:
            conflicts.setdefault(rank_name, [taxonomy[rank_name]]).append(value)
        taxonomy[rank_name] = value if value is not None else None
        if value is None:
            self.logger.debug(f"Empty value for rank '{rank_name}' in taxonomy string")

    def save_to_file(self, file_path: Path, names: Optional[List[str]] = None) -> int:
        """
        保存分类映射到 TSV 文件

        Args:
            file_path: 输出文件路径
            names: 要保存的名称列表 (None = 所有)

        Returns:
            实际写出的条目数；``0`` 表示**没有写出任何文件**（调用方可据此判断
            产物是否真的存在）。

        Note:
            指定 ``names`` 时与 ``names=None`` 走同一优先序的数据源：先查已加载的
            外部分类表，再查解析缓存（C-44——旧实现只查 ``_cache``，于是
            "刚 load_from_file 进来的表"一个都写不出去，只留一条 warning 就产出
            空产物）。查不到的名字逐个披露并计数。
        """
        all_entries: Dict[str, TaxonomyMap] = {}
        if names is None:
            for name, taxonomy in list(self._file_loaded_taxonomy.items()) + list(
                self._cache.items()
            ):
                if name not in all_entries:
                    all_entries[name] = taxonomy
        else:
            missing: List[str] = []
            for name in names:
                # 这里的变量不能和上面的 ``taxonomy`` 同名：循环变量在两个分支之间
                # 共享一个类型，而本分支拿到的是 ``Optional[TaxonomyMap]``。
                entry = self._file_loaded_taxonomy.get(name)
                if entry is None:
                    entry = self._cache.get(name)
                if entry is None:
                    missing.append(name)
                    continue
                all_entries[name] = entry
            if missing:
                self.logger.warning(
                    f"{len(missing)} of {len(names)} requested name(s) have no "
                    f"taxonomy entry and were skipped: "
                    f"{missing[:10]}{'...' if len(missing) > 10 else ''}"
                )

        if not all_entries:
            self.logger.warning(
                f"No taxonomy entries to save"
                f"{'' if names is None else f' (asked for {len(names)} name(s))'}; "
                f"{file_path} was NOT written."
            )
            return 0

        file_path.parent.mkdir(parents=True, exist_ok=True)

        ranks = list(self._levels.values())

        with safe_writer(file_path, newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, delimiter="\t", fieldnames=["name"] + ranks)
            writer.writeheader()

            for name, taxonomy in all_entries.items():
                row = {"name": name}
                for rank in ranks:
                    value = taxonomy.get(rank, "")
                    row[rank] = value if value is not None else ""
                writer.writerow(row)

        self.logger.info(f"Saved {len(all_entries)} taxonomy entries to {file_path}")
        return len(all_entries)

    def get_taxonomy_source(self, name: str) -> Optional[str]:
        """
        获取分类信息来源

        Returns:
            'file', 'gtdb', 'gtdb_underscore', 'ncbi', 或 None
        """
        return self._taxonomy_source.get(name)

    @property
    def has_external_file(self) -> bool:
        """是否加载了外部分类文件"""
        return len(self._file_loaded_taxonomy) > 0

    @property
    def external_file_path(self) -> Optional[Path]:
        """获取外部分类文件路径"""
        return self._loaded_file_path

    @property
    def supported_levels(self) -> Dict[str, str]:
        """获取支持的分类级别映射"""
        return self._levels.copy()

    def detect_format(self, name: str) -> Optional[str]:
        """
        检测分类学字符串的格式类型

        Args:
            name: 叶节点名称或分类学字符串

        Returns:
            格式类型: 'format_a', 'format_b', 'ncbi', 或 None
        """
        if not name:
            return None

        if "__" in name and ";" in name:
            return "format_b"

        if re.search(r"_[a-z]_", name):
            return "format_a"

        if ";" in name and not any(c in name for c in ["_", "__"]):
            return "ncbi"

        return None

    def merge_taxonomy(self, *taxonomy_dicts: Optional[TaxonomyMap]) -> TaxonomyMap:
        """
        合并多个分类学字典

        用于统一来自不同格式的分类学信息。
        后面的字典会覆盖前面的字典中的相同键。

        Args:
            *taxonomy_dicts: 多个分类学字典

        Returns:
            合并后的分类学字典
        """
        merged: TaxonomyMap = {}

        for taxonomy in taxonomy_dicts:
            if taxonomy is None:
                continue

            for rank, value in taxonomy.items():
                if value is not None and value:
                    merged[rank] = value
                elif rank not in merged:
                    merged[rank] = None

        return merged

    def format_taxonomy_string(
        self, taxonomy: TaxonomyMap, format_type: str = "format_b"
    ) -> str:
        """
        将分类学字典格式化为字符串

        Args:
            taxonomy: 分类学字典
            format_type: 输出格式 ('format_a', 'format_b', 'ncbi')

        Returns:
            格式化的分类学字符串
        """
        if not taxonomy:
            return ""

        if format_type == "format_b":
            parts = []
            for rank_name, value in taxonomy.items():
                prefix = self._rank_to_prefix.get(rank_name)
                if prefix:
                    value_str = value if value is not None else ""
                    parts.append(f"{prefix}__{value_str}")
            return ";".join(parts)

        elif format_type == "format_a":
            parts = []
            for rank_name, value in taxonomy.items():
                prefix = self._rank_to_prefix.get(rank_name)
                if prefix:
                    value_str = value if value is not None else ""
                    parts.append(f"_{prefix}_{value_str}")
            return "".join(parts)

        elif format_type == "ncbi":
            parts = []
            for rank_name in [
                "domain",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ]:
                value = taxonomy.get(rank_name)
                if value:
                    parts.append(value)
            return ";".join(parts)

        return ""
