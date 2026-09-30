"""
NameMappingManager - 名称映射管理器

管理原始分类学名称与简化名称之间的双向映射

映射层的核心契约（审阅报告 B-21）：**短名唯一性是双向映射成立的前提**。
任一时刻都必须满足 ``len(get_reverse_mapping_dict()) == len(get_mapping_dict())``，
否则反查（``_restore_leaf_names`` 一类）会把两个分类单元还原成同一个全名——
一个从输出树里消失、另一个出现两次，而输出树看上去完全正常。
因此本模块：

- 生成短名时一律经 ``resolve_short_name`` 做冲突检测，冲突时追加稳定的
  计数后缀（``_2`` / ``_3`` …）并在长度上限内截断，绝不静默覆盖；
- 反查表接口默认以 ``strict=True`` 自检完整性，一旦发现"有分类单元在反查时
  不可见"直接 ``NameCollisionError`` 上抛（``strict=False`` 时降级为 warning，
  但仍然不丢条目）；
- 从文件加载映射时同样拒绝同短名 / 同全名的静默覆盖。
"""

import csv
import re
import threading
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Iterable, List, Optional, OrderedDict

from ..infrastructure.logging import get_logger
from ..infrastructure.safe_io import safe_writer


class NameShortenMode(Enum):
    """名称缩短模式"""

    SERIAL = "serial"  # 串行编号（≤10字符，MCMCTree严格PHYLIP）
    ACCESSION = "accession"  # 登录号（≤30字符，PAML标准PHYLIP）
    FULL = "full"  # 完整名称（无限制）


# 各模式短名的长度上限（None = 不限制）。消歧后缀必须在这个预算之内生成，
# 否则"加了后缀就不超限"的短名会在下游（PHYLIP 定长列）再次被截断成同名。
SHORT_NAME_LENGTH_LIMITS: Dict[NameShortenMode, Optional[int]] = {
    NameShortenMode.SERIAL: 10,
    NameShortenMode.ACCESSION: 30,
    NameShortenMode.FULL: None,
}


class NameCollisionError(ValueError):
    """短名冲突无法消歧，或反查表会静默丢失分类单元（B-21）。"""


def _fit_with_suffix(base: str, suffix: str, max_length: Optional[int]) -> str:
    """把 ``base`` + ``suffix`` 压进 ``max_length``（None 表示不限制）。"""
    if max_length is None:
        return f"{base}{suffix}"
    keep = max_length - len(suffix)
    if keep <= 0:
        return suffix[:max_length]
    return f"{base[:keep]}{suffix}"


def resolve_short_name(
    base: str,
    claimed: Dict[str, str],
    owner: str,
    max_length: Optional[int] = None,
    max_attempts: int = 100000,
) -> str:
    """为 ``owner`` 申请一个未被占用的短名（B-21 的唯一性检查）。

    Args:
        base: 期望的短名（如登录号、串行号、清洗后的全名）
        claimed: 已占用的 ``短名 -> 归属原始名`` 表
        owner: 本次申请的原始名（同名的重复申请是幂等的）
        max_length: 短名长度上限，None 表示不限制
        max_attempts: 消歧后缀的最大尝试次数，超过即判定无法消歧

    Returns:
        唯一短名：``base``（在 ``max_length`` 内）未被子占用时原样返回；
        否则返回 ``base_2`` / ``base_3`` … 形式的消歧结果。

    Raises:
        NameCollisionError: 后缀空间耗尽仍无法找到唯一短名
    """
    plain = base if max_length is None else base[:max_length]
    if not plain:
        raise NameCollisionError(
            f"无法为 {owner!r} 生成短名：期望短名 {base!r} 为空"
            f"（或在长度上限 {max_length} 下被截为空串）。"
        )
    candidate = plain
    holder = 1
    tried = {plain}
    while True:
        current = claimed.get(candidate)
        if current is None or current == owner:
            return candidate
        holder += 1
        if holder > max_attempts:
            raise NameCollisionError(
                f"短名冲突无法消歧: {base!r}（归属 {owner!r}）在 {max_attempts} 次"
                f"追加后缀后仍与既有短名相撞，请缩短名称或改用其它缩短模式。"
            )
        candidate = _fit_with_suffix(base, f"_{holder}", max_length)
        if candidate in tried:
            # 长度预算太小，后缀本身已被截成同一个串：继续下去只是原地打转
            raise NameCollisionError(
                f"短名冲突无法消歧: 长度上限 {max_length} 下 {base!r} 的消歧后缀"
                f"已用完（{candidate!r} 重复），请放宽上限或改用其它缩短模式。"
            )
        tried.add(candidate)


# sanitize_name 的字符清洗规则（同一份正则被 sanitize_name 与批量入口共用，
# 避免"两份拷贝改一份"导致的树/约束文件名不一致）。
_SANITIZE_PATTERN = re.compile(r"[\s\(\)\[\]\{\}\|\'\"\;:,]+")
_COLLAPSE_PATTERN = re.compile(r"_+")


def sanitize_name_basic(name: str, max_length: Optional[int] = 30) -> str:
    """清洗单个名称中的特殊字符并按需截断（不消歧，见 ``sanitize_names_uniquely``）。"""
    sanitized = _SANITIZE_PATTERN.sub("_", name)
    sanitized = _COLLAPSE_PATTERN.sub("_", sanitized)
    sanitized = sanitized.strip("_")
    if max_length is not None and len(sanitized) > max_length:
        sanitized = sanitized[:max_length]
    return sanitized


def sanitize_names_uniquely(
    names: Iterable[str],
    max_length: int = 30,
) -> Dict[str, str]:
    """批量清洗名称并保证结果两两不同（B-21 / C-30 的共用注册入口）。

    ``sanitize_name`` 是纯函数：截断到 30 字符时两个不同全名可以塌成同一个
    短名（审阅报告 C-30 的实测例）。各适配器（r8s / treePL / PATHd8）此前各自
    实现 ``_sanitize_name`` 而不过 ``NameMappingManager``，于是连"两个表对不上"
    的机会都没有。本函数把"清洗 + 按输入顺序追加 ``_2`` / ``_3`` … 后缀"合并为
    一个确定性入口，供适配器直接替换裸 ``[:30]``。

    Args:
        names: 原始名称（可迭代，允许重复；重复者共享同一短名）
        max_length: 短名长度上限（含消歧后缀）

    Returns:
        ``原始名 -> 唯一短名`` 字典，按输入顺序生成，因此同名映射是确定性的
    """
    claimed: Dict[str, str] = {}
    mapping: Dict[str, str] = {}
    for index, name in enumerate(names, start=1):
        if name in mapping:
            continue
        # 先完整清洗（不截断），再由 resolve_short_name 在长度预算内截断 + 消歧，
        # 保证后缀一定落在 max_length 之内。整串都是被替换字符时用序号兜底，
        # 避免下游拿到空叶名。
        base = sanitize_name_basic(name, max_length=None) or f"seq{index:05d}"
        short = resolve_short_name(base, claimed, owner=name, max_length=max_length)
        claimed[short] = name
        mapping[name] = short
    return mapping


@dataclass
class NameMappingEntry:
    """名称映射条目"""

    original: str  # 原始名称
    short: str  # 简化名称
    accession: str  # 登录号
    taxonomy: Dict  # 分类学信息
    method: NameShortenMode  # 重命名策略


class NameMappingManager:
    """
    名称映射管理器

    管理原始分类学名称与简化名称之间的双向映射
    线程安全
    """

    def __init__(self) -> None:
        self._entries: OrderedDict[str, NameMappingEntry] = OrderedDict()
        self._short_to_original: Dict[str, str] = {}
        # 登录号 -> 原始名。共享同一登录号的多个序列（GTDB 风格的同一 assembly
        # 拆成多条分片）在单值字典里会互相遮蔽，因此主存储是"多值表"；
        # ``_accession_to_original`` 仅作先入者的快捷视图保留（B-21）。
        self._accession_to_original: Dict[str, str] = {}
        self._accession_to_originals: Dict[str, List[str]] = {}
        self._counter = 0
        # 线程锁，确保线程安全
        self._lock = threading.Lock()
        self.logger = get_logger()

        # 实例级解析缓存：分类学解析与登录号提取均为名称的纯函数，对同一名称
        # 仅解析一次后缓存（#16 性能/可维护性收敛），避免重复 IO 与正则解析。
        self._taxonomy_cache: Dict[str, Dict[str, str]] = {}
        self._accession_cache: Dict[str, str] = {}

    def add_entry(
        self,
        original: str,
        method: NameShortenMode = NameShortenMode.ACCESSION,
        on_method_mismatch: str = "warn",
    ) -> str:
        """
        添加名称映射条目（线程安全、短名唯一）

        Args:
            original: 原始名称
            method: 重命名策略
            on_method_mismatch: 同一 original 以不同 method 重复注册时的处理
                （C-29）：``'warn'``（默认）记 WARNING 后沿用已生成的短名；
                ``'raise'`` 抛 ``ValueError``；``'regenerate'`` 按新策略重新生成
                短名并做冲突消歧（旧短名一并让位，反查表保持一一对应）。

        Returns:
            简化后的名称（保证与本管理器内其它条目的短名互不相同）

        Raises:
            ValueError: 如果 original 为空或 None，或 ``on_method_mismatch='raise'``
                且本次 method 与已注册条目不一致
            TypeError: 如果 original 不是字符串
            NameCollisionError: 短名冲突且无法消歧
        """
        # 空值验证
        if original is None:
            raise TypeError("original name cannot be None")
        if not isinstance(original, str):
            raise TypeError(
                f"original name must be a string, got {type(original).__name__}"
            )
        if not original.strip():
            raise ValueError("original name cannot be empty or whitespace only")
        if on_method_mismatch not in ("warn", "raise", "regenerate"):
            raise ValueError(
                f"on_method_mismatch must be 'warn', 'raise' or 'regenerate', "
                f"got {on_method_mismatch!r}"
            )

        with self._lock:
            if original in self._entries:
                # C-29：调用方本次要求的策略与已注册条目不一致时不再静默丢弃
                self._handle_method_mismatch(original, method, on_method_mismatch)
                if original in self._entries:
                    return self._entries[original].short

            # 解析登录号
            accession = self._extract_accession(original)

            # 解析分类学信息
            taxonomy = self._parse_taxonomy(original)

            # 生成简化名称
            if method == NameShortenMode.SERIAL:
                base = self._generate_serial_name()
            elif method == NameShortenMode.ACCESSION:
                base = accession if accession else self._generate_serial_name()
            else:  # FULL
                base = original

            # B-21：短名必须唯一——登录号不是序列的唯一标识（同一 GCA_… 的多个
            # 分片在 GTDB 风格数据里很常见），直接写入反查表会让前一个分类单元
            # 无声消失。冲突时追加稳定计数后缀，长度仍受本模式上限约束。
            short = self._claim_short_name(base, original, method)

            entry = NameMappingEntry(
                original=original,
                short=short,
                accession=accession,
                taxonomy=taxonomy,
                method=method,
            )

            self._entries[original] = entry
            if accession:
                bucket = self._accession_to_originals.setdefault(accession, [])
                if original not in bucket:
                    bucket.append(original)
                # 先入者优先（与 taxonomy_parser.save_to_file / merge_taxonomy 一致）
                self._accession_to_original.setdefault(accession, original)

            return short

    def _handle_method_mismatch(
        self, original: str, method: NameShortenMode, on_method_mismatch: str
    ) -> None:
        """处理"同一 original、不同 method"的重复注册（C-29）。

        调用本方法时已持有 ``self._lock``。``'regenerate'`` 分支会就地重建短名，
        因此调用方需要在返回后重新检查 ``original`` 是否仍在 ``_entries`` 中。
        """
        existing = self._entries.get(original)
        if existing is None or existing.method == method:
            return

        if on_method_mismatch == "raise":
            raise ValueError(
                f"名称 {original!r} 已按 method={existing.method.value} 注册为短名 "
                f"{existing.short!r}，本次请求 method={method.value} 被拒绝"
                f"（on_method_mismatch='raise'）。"
            )

        if on_method_mismatch == "warn":
            self.logger.warning(
                f"add_entry 策略冲突: {original!r} 已按 "
                f"method={existing.method.value} 注册为 {existing.short!r}，"
                f"本次 method={method.value} 被忽略（沿用已生成的短名）。"
                f"如需按新策略重建请传 on_method_mismatch='regenerate'。"
            )
            return

        # 'regenerate'：撤销旧条目并按新策略重建（旧短名让位，避免反查表残留）
        old_short = existing.short
        if self._short_to_original.get(old_short) == original:
            del self._short_to_original[old_short]
        del self._entries[original]
        self.logger.debug(
            f"add_entry 按 method={method.value} 重新生成 {original!r} 的短名"
            f"（原 {old_short!r}）"
        )

    def _claim_short_name(
        self, base: str, original: str, method: NameShortenMode
    ) -> str:
        """申请唯一短名并登记反查表（调用方须持有 ``self._lock``）。

        与既有短名冲突时不覆盖、不改数据：追加 ``_2`` / ``_3`` … 后缀并记 WARNING，
        消歧空间不足则由 ``resolve_short_name`` 抛 ``NameCollisionError``。
        """
        limit = SHORT_NAME_LENGTH_LIMITS.get(method)
        plain = base if limit is None else base[:limit]
        previous_holder = self._short_to_original.get(plain)
        short = resolve_short_name(
            plain, self._short_to_original, owner=original, max_length=limit
        )
        if short != plain:
            self.logger.warning(
                f"短名冲突已消歧: {plain!r} 已被 "
                f"{previous_holder!r} 占用，{original!r} 改用唯一短名 "
                f"{short!r}（B-21：不消歧会让反查表静默丢掉一个分类单元）。"
            )
        self._short_to_original[short] = original
        return short

    def _extract_accession(self, name: str) -> str:
        """从名称中提取登录号（结果按名称实例级缓存，避免重复解析）"""
        cached = self._accession_cache.get(name)
        if cached is not None:
            return cached
        # GTDB 格式: GB_GCA_000252485.1_d_Bacteria...
        match = re.search(r"([A-Z]{2}_[A-Z]{2,3}_\d+\.\d+)", name)
        if match:
            accession = match.group(1).replace(".", "_")
        else:
            # 其他格式尝试提取
            match = re.search(r"(GC[AF]_\d+\.\d+)", name)
            accession = match.group(1).replace(".", "_") if match else ""
        self._accession_cache[name] = accession
        return accession

    def _parse_taxonomy(self, name: str) -> Dict[str, str]:
        """从名称中解析分类学信息（结果按名称实例级缓存，避免重复解析）"""
        cached = self._taxonomy_cache.get(name)
        if cached is not None:
            return cached
        taxonomy = {}

        # GTDB 下划线分隔格式
        patterns = [
            (r"_d_(?P<domain>[^_]+)", "domain"),
            (r"_p_(?P<phylum>[^_]+)", "phylum"),
            (r"_c_(?P<class>[^_]+)", "class"),
            (r"_o_(?P<order>[^_]+)", "order"),
            (r"_f_(?P<family>[^_]+)", "family"),
            (r"_g_(?P<genus>[^_]+)", "genus"),
            (r"_s_(?P<species>[^_]+)", "species"),
        ]

        for pattern, rank in patterns:
            match = re.search(pattern, name)
            if match:
                taxonomy[rank] = match.group(1)

        self._taxonomy_cache[name] = taxonomy
        return taxonomy

    def _generate_serial_name(self) -> str:
        """生成串行编号名称"""
        self._counter += 1
        return f"S{self._counter:05d}"

    def get_short_name(self, original: str) -> Optional[str]:
        """获取简化名称（线程安全）"""
        with self._lock:
            entry = self._entries.get(original)
            return entry.short if entry else None

    def get_original_name(self, short: str) -> Optional[str]:
        """从简化名称获取原始名称（线程安全）"""
        with self._lock:
            return self._short_to_original.get(short)

    def get_accession(self, original: str) -> Optional[str]:
        """获取登录号（线程安全）"""
        with self._lock:
            entry = self._entries.get(original)
            return entry.accession if entry else None

    def get_originals_by_accession(self, accession: str) -> List[str]:
        """获取共享同一登录号的**全部**原始名（线程安全）。

        登录号不是序列的唯一标识（B-21），因此这里必须是多值表：任何按登录号
        反查"是哪个分类单元"的代码都应当先问有几个候选，而不是取第一个。
        """
        with self._lock:
            return list(self._accession_to_originals.get(accession, ()))

    def get_taxonomy(self, original: str) -> Optional[Dict[str, str]]:
        """获取分类学信息（线程安全）"""
        with self._lock:
            entry = self._entries.get(original)
            return entry.taxonomy if entry else None

    def get_mapping_dict(self) -> Dict[str, str]:
        """获取原始名称到简化名称的映射字典（线程安全）"""
        with self._lock:
            return {entry.original: entry.short for entry in self._entries.values()}

    def find_duplicate_short_names(self) -> Dict[str, List[str]]:
        """列出被两个以上原始名共用的短名（线程安全，正常应恒为空）。

        返回 ``{短名: [占用它的原始名, …]}``；非空即意味着反查表会丢掉分类单元。
        """
        with self._lock:
            return self._find_duplicate_short_names()

    def _find_duplicate_short_names(self) -> Dict[str, List[str]]:
        """``find_duplicate_short_names`` 的无锁实现（调用方须持有 ``self._lock``）。"""
        groups: Dict[str, List[str]] = {}
        for entry in self._entries.values():
            groups.setdefault(entry.short, []).append(entry.original)
        return {s: names for s, names in groups.items() if len(names) > 1}

    def check_mapping_integrity(self) -> List[str]:
        """自检双向映射是否一一对应，返回问题描述列表（线程安全，无副作用）。

        用于任何一次批量映射之后（审阅报告 B-21 的修复建议：
        ``len(set(mapping.values())) == len(mapping)`` 断言）。
        """
        with self._lock:
            return self._check_mapping_integrity()

    def _check_mapping_integrity(self) -> List[str]:
        """``check_mapping_integrity`` 的无锁实现（调用方须持有 ``self._lock``）。"""
        problems: List[str] = []

        for short, names in self._find_duplicate_short_names().items():
            problems.append(f"短名 {short!r} 被 {len(names)} 个原始名共用: {names}")

        for entry in self._entries.values():
            holder = self._short_to_original.get(entry.short)
            if holder != entry.original:
                problems.append(
                    f"反查表不可达: 短名 {entry.short!r} 现指向 {holder!r}，"
                    f"而非其归属 {entry.original!r}"
                )

        if len(self._short_to_original) != len(self._entries):
            problems.append(
                f"条目数不一致: 正向 {len(self._entries)} 条 vs 反向 "
                f"{len(self._short_to_original)} 条"
            )

        return problems

    def assert_mapping_integrity(self) -> None:
        """双向映射不完整时直接 ``NameCollisionError`` 上抛（线程安全）。"""
        problems = self.check_mapping_integrity()
        if problems:
            raise NameCollisionError(
                "名称映射完整性检查失败，反查会静默丢失分类单元"
                "（输出树里会有一个分类单元消失、另一个出现两次）：\n  - "
                + "\n  - ".join(problems)
            )

    def get_reverse_mapping_dict(self, strict: bool = True) -> Dict[str, str]:
        """获取简化名称到原始名称的映射字典（线程安全）。

        Args:
            strict: True（默认）时先自检完整性，发现任何分类单元在反查表中不可见
                即抛 ``NameCollisionError``——静默丢掉一个 taxon 比失败危险得多
                （B-21）。False 时降级为 WARNING 并返回，供只做诊断的调用方使用。
        """
        with self._lock:
            problems = self._check_mapping_integrity()
            if problems:
                message = (
                    "名称映射完整性检查失败，反查会静默丢失分类单元："
                    + "; ".join(problems)
                )
                if strict:
                    raise NameCollisionError(message)
                self.logger.warning(message)
            return dict(self._short_to_original)

    def save_to_file(self, file_path: Path) -> None:
        """保存映射表到 TSV 文件（线程安全）"""
        with self._lock:
            problems = self._check_mapping_integrity()
            if problems:
                raise NameCollisionError(
                    f"拒绝写出名称映射 {file_path}：短名不唯一会使下游反查静默丢失"
                    f"分类单元（B-21）。问题：\n  - " + "\n  - ".join(problems)
                )

            file_path.parent.mkdir(parents=True, exist_ok=True)

            with safe_writer(file_path, newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(
                    ["original", "short", "accession", "taxonomy", "method"]
                )

                for entry in self._entries.values():
                    taxonomy_str = ";".join(
                        f"{k}={v}" for k, v in entry.taxonomy.items()
                    )
                    writer.writerow(
                        [
                            entry.original,
                            entry.short,
                            entry.accession,
                            taxonomy_str,
                            entry.method.value,
                        ]
                    )

    def load_from_file(self, file_path: Path) -> None:
        """从 TSV 文件加载映射表（线程安全，拒绝会丢失分类单元的映射）。

        文件里的 ``short`` 列若在两行之间重复，反查表就会静默丢掉一个 taxon
        （B-21 的持久化形态）——这里直接 ``NameCollisionError`` 上抛，而不是
        "后一行覆盖前一行"。校验失败时**原有映射保持不变**（先在局部字典里
        构建，成功后才提交），避免"加载坏文件 → 映射被清空"的连带失效。

        Raises:
            ValueError: 表头缺失 / 必需列缺失 / 未知 method 值
            NameCollisionError: 同一 short 被两行占用，或同一 original 出现两行
                且短名不一致
        """
        with self._lock:
            entries: "OrderedDict[str, NameMappingEntry]" = OrderedDict()
            short_to_original: Dict[str, str] = {}
            accession_to_original: Dict[str, str] = {}
            accession_to_originals: Dict[str, List[str]] = {}

            with open(file_path, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")

                # 验证必需列
                required = {"original", "short", "method"}
                if reader.fieldnames is None:
                    raise ValueError(f"Mapping file {file_path} has no header row")
                missing = required - set(reader.fieldnames)
                if missing:
                    raise ValueError(
                        f"Mapping file {file_path} is missing required columns: {sorted(missing)}"
                    )

                for line_no, row in enumerate(reader, start=2):
                    taxonomy = {}
                    taxonomy_raw = (row.get("taxonomy") or "").strip()
                    if taxonomy_raw:
                        for part in taxonomy_raw.split(";"):
                            if "=" in part:
                                k, v = part.split("=", 1)
                                taxonomy[k] = v

                    accession = (row.get("accession") or "").strip()
                    original = row["original"]
                    short = row["short"]

                    previous_original = short_to_original.get(short)
                    if previous_original is not None and previous_original != original:
                        raise NameCollisionError(
                            f"映射文件 {file_path} 第 {line_no} 行：短名 {short!r} 已被 "
                            f"{previous_original!r} 使用，本行却属于 {original!r}。"
                            f"继续加载会让反查表只剩一个归属、另一个分类单元在还原时"
                            f"消失（B-21）。请重新生成该映射文件。"
                        )
                    if original in entries and entries[original].short != short:
                        raise NameCollisionError(
                            f"映射文件 {file_path} 第 {line_no} 行：{original!r} 出现多次"
                            f"且短名不一致（{entries[original].short!r} vs {short!r}），"
                            f"无法确定唯一的还原目标。"
                        )

                    entry = NameMappingEntry(
                        original=original,
                        short=short,
                        accession=accession,
                        taxonomy=taxonomy,
                        method=NameShortenMode(row["method"]),
                    )

                    entries[original] = entry
                    short_to_original[short] = original
                    if accession:
                        accession_to_original.setdefault(accession, original)
                        bucket = accession_to_originals.setdefault(accession, [])
                        if original not in bucket:
                            bucket.append(original)

            # 通过全部校验后才提交
            self._entries = entries
            self._short_to_original = short_to_original
            self._accession_to_original = accession_to_original
            self._accession_to_originals = accession_to_originals

            # 从已有短名中恢复计数器，避免后续生成重复短名。
            # 同时识别消歧后缀形态（S00007_2），否则计数器会退回已用编号，
            # 新生成的短名虽仍会被唯一性检查救下，但编号会出现回绕。
            nums = [
                int(m.group(1))
                for s in self._short_to_original
                if (m := re.fullmatch(r"S(\d{5})(?:_\d+)?", s))
            ]
            self._counter = max(nums, default=0)

    def get_entries_by_taxonomy(self, rank: str, value: str) -> List[str]:
        """
        根据分类学信息获取原始名称列表（线程安全）

        Args:
            rank: 分类级别（domain, phylum, class, order, family, genus）
            value: 分类值

        Returns:
            匹配的原始名称列表
        """
        with self._lock:
            result = []
            for entry in self._entries.values():
                if entry.taxonomy.get(rank) == value:
                    result.append(entry.original)
            return result

    @staticmethod
    def sanitize_name(name: str, max_length: int = 30) -> str:
        """
        清洗序列名称中的特殊字符

        替换以下字符为空格然后压缩:
        - 空格、制表符
        - 括号 ()[]{}
        - 竖线 |
        - 单引号、双引号
        - 分号、冒号
        - 逗号

        Args:
            name: 原始名称
            max_length: 最大长度限制

        Returns:
            清洗后的安全名称

        Note:
            本方法是**纯函数**，只清洗单个名称，因此**不保证唯一**：截断到
            ``max_length`` 时两个不同全名会塌成同一个短名（C-30，实测
            ``Escherichia_coli_strain_ABCDEFGHIJ_more_sequence_data`` 与
            ``…_other_sequence_data`` 同为 30 字符前缀）。批量注册请使用
            ``sanitize_names_uniquely``（清洗 + 冲突消歧）或
            ``NameMappingManager.add_entry``，不要各自 ``[:30]``。
        """
        return sanitize_name_basic(name, max_length=max_length)
