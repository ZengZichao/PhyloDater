"""
CalibrationLoader - 校准配置加载器

从 YAML 文件加载校准点配置

本模块同时是**校准 YAML 内容契约的唯一权威**（审阅项 B-24 / C-38）：
``CONSTRAINT_SPECS`` 描述八类合法约束的必填/可选键与取值域，
:func:`validate_calibration_payload` 逐条目做内容级校验，
:class:`~phylodater.services.input_validator.InputValidator` 的
``validate_calibration_file`` 与本加载器共用这一份表，避免两处漂移。
"""

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

from ..core.exceptions import CalibrationError
from ..infrastructure import get_logger
from ..infrastructure.safe_io import safe_writer
from ..models import (
    CalibrationPoint,
    FixedAgeConstraint,
    FossilConfidence,
    FossilMetadata,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    NodeDefinition,
    SkewNormalConstraint,
    SkewTConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)
from ..models.constraints import MAX_GEOLOGICAL_AGE_MA

# --------------------------------------------------------------------------- #
# 校准 YAML 的内容契约（B-24）
# --------------------------------------------------------------------------- #

# 旧/别名约束类型映射（向后兼容）；键为 YAML 里写的名字，值为规范类型名
CONSTRAINT_TYPE_ALIASES: Dict[str, str] = {
    "lower_bound": "soft_lower",
    "upper_bound": "maximum",
    "point": "fixed",
    "range": "uniform",
}

# 单个约束块的键规格：
#   required     : 必须存在的键（缺失即条目无效，不再"猜"成 fixed）
#   optional     : 允许出现的键（其余键报 warning，因为它们通常意味着拼写错误）
#   age_keys     : 以 Ma 计的年龄键，必须是有限、0 <= v <= 4600
#   prob_keys    : 概率键，必须 0 < v < 1
#   positive_keys: 必须为正的有限数值（尺度/速率/形状类）
#   note         : 面向用户的语义提示（同一个 YAML 键在不同类型下语义不同时必须说清）
CONSTRAINT_SPECS: Dict[str, Dict[str, Any]] = {
    "fixed": {
        "required": ("age",),
        "optional": (),
        "age_keys": ("age",),
        "prob_keys": (),
        "positive_keys": (),
    },
    "uniform": {
        "required": ("min", "max"),
        "optional": ("p_lower", "p_upper"),
        "age_keys": ("min", "max"),
        "prob_keys": ("p_lower", "p_upper"),
        "positive_keys": (),
        "note": "uniform 是近似硬界；'p' 对它没有意义（要用软界请写 type: soft_bounds）",
    },
    "soft_lower": {
        "required": ("min",),
        "optional": ("p", "c", "tail_prob"),
        "age_keys": ("min",),
        "prob_keys": ("tail_prob",),
        "positive_keys": ("c",),
        # 审阅项 B-2：YAML 的 p/c 就是 PAML 手册里 L(tL, p, c) 的 p 与 c，
        # 即**相对偏移比例**与 **Cauchy 尺度**，不是尾概率。
        "note": "'p' 是 MCMCTree L() 的相对偏移 P（先验中心 = min*(1+P)），"
        "不是概率；'c' 是 Cauchy 尺度；真正的左尾概率是 'tail_prob'（L() 槽 4）",
    },
    "maximum": {
        "required": ("max",),
        "optional": ("p",),
        "age_keys": ("max",),
        "prob_keys": ("p",),
        "positive_keys": (),
    },
    "soft_bounds": {
        "required": ("min", "max"),
        "optional": ("p", "p_lower", "p_upper"),
        "age_keys": ("min", "max"),
        "prob_keys": ("p", "p_lower", "p_upper"),
        "positive_keys": (),
    },
    "gamma": {
        "required": ("alpha", "beta"),
        "optional": ("offset", "scale"),
        "age_keys": ("offset",),
        "prob_keys": (),
        "positive_keys": ("alpha", "beta", "scale"),
    },
    "skew_normal": {
        "required": ("location", "scale", "shape"),
        "optional": (),
        "age_keys": ("location",),
        "prob_keys": (),
        "positive_keys": ("scale",),
    },
    "skew_t": {
        "required": ("location", "scale", "shape", "df"),
        "optional": (),
        "age_keys": ("location",),
        "prob_keys": (),
        "positive_keys": ("scale", "df"),
    },
}

# ``df`` 除了为正还必须 > 2（t 分布方差有限的前提）
_MIN_DF = 2.0


def constraint_type_names() -> List[str]:
    """全部合法约束类型名（含别名），用于错误提示文案。"""
    return list(CONSTRAINT_SPECS) + [
        f"{a} (->{t})" for a, t in CONSTRAINT_TYPE_ALIASES.items()
    ]


def format_constraint_help() -> str:
    """一行式的约束类型/字段提示（CLI 与库层共用，避免三处各写一份）。"""
    fields = sorted(
        {
            k
            for spec in CONSTRAINT_SPECS.values()
            for k in spec["required"] + spec["optional"]
        }
    )
    return (
        "支持的约束类型 ("
        + "/".join(CONSTRAINT_SPECS)
        + "，别名 "
        + "/".join(CONSTRAINT_TYPE_ALIASES)
        + ") 以及合法字段 ("
        + "/".join(fields)
        + ")"
    )


def _canonical_type_name(raw: Any) -> Optional[str]:
    """把 YAML 写的类型名规范成 ``CONSTRAINT_SPECS`` 的键（别名可被识别）。"""
    if not isinstance(raw, str):
        return None
    name = raw.strip().lower()
    if name in CONSTRAINT_SPECS:
        return name
    return CONSTRAINT_TYPE_ALIASES.get(name)


def _as_finite_float(
    value: Any, key: str, where: str
) -> Tuple[Optional[float], Optional[str]]:
    """把 YAML 值转成**有限**浮点数（B-19：拒绝 ``.nan`` / ``.inf`` / ``nan`` 字符串）。

    Returns:
        (数值, None) 或 (None, 错误信息)
    """
    if isinstance(value, bool):
        return None, f"{where}: 键 '{key}' 必须是数值，收到布尔值 {value!r}"
    if not isinstance(value, (int, float, str)):
        return None, (
            f"{where}: 键 '{key}' 必须是数值，收到 {type(value).__name__} {value!r}"
        )
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None, f"{where}: 键 '{key}' 无法解析为数值: {value!r}"
    if math.isnan(number):
        return None, (
            f"{where}: 键 '{key}' 是 NaN（YAML 的 .nan 字面量）。NaN 与任何数的比较"
            f"都为 False，会静默绕过所有上下界检查（审阅项 B-19），因此被拒绝"
        )
    if math.isinf(number):
        return None, f"{where}: 键 '{key}' 是无穷大（{value!r}），必须是有限数值"
    return number, None


def validate_constraint_block(data: Any, where: str = "constraint") -> List[str]:
    """内容级校验单个 ``constraint`` 块，返回错误信息列表（空 = 通过）。"""
    errors: List[str] = []

    if not isinstance(data, dict):
        return [
            f"{where}: 必须是一个映射（含 type/age/min/max 等键），"
            f"实际收到 {type(data).__name__}: {data!r}"
        ]

    if "type" not in data:
        # 审阅项 C-38：缺 type 时旧实现静默按 fixed 处理，于是真正的问题
        # （type 放错层级 / 忘写）被翻译成"缺 age 键"这种指错方向的报错。
        return [
            f"{where}: 缺少 'type' 键。已知键为 {sorted(map(str, data.keys()))}。"
            f"缺 type 时不再猜测（旧版本会按 'fixed' 处理），{format_constraint_help()}"
        ]

    raw_type = data["type"]
    canonical = _canonical_type_name(raw_type)
    if canonical is None:
        return [
            f"{where}: 未知的约束类型 'type: {raw_type}'。"
            f"{format_constraint_help()}"
        ]

    spec = CONSTRAINT_SPECS[canonical]
    alias_note = ""
    if (
        isinstance(raw_type, str)
        and raw_type.strip().lower() in CONSTRAINT_TYPE_ALIASES
    ):
        alias_note = f"（由旧类型名 '{raw_type}' 翻译而来）"

    for key in spec["required"]:
        if key not in data:
            errors.append(
                f"{where}: 类型 '{canonical}'{alias_note} 缺少必填键 '{key}' "
                f"(已知键 {sorted(map(str, data.keys()))})"
            )

    allowed = set(spec["required"]) | set(spec["optional"])
    for key in data:
        if key == "type":
            continue
        if key not in allowed:
            errors.append(
                f"{where}: 类型 '{canonical}'{alias_note} 不认识键 '{key}'"
                f"（可用键: {sorted(allowed)}）。拼错的键会被忽略，"
                f"导致你以为生效的参数实际没有生效"
            )

    for key in spec["age_keys"]:
        if key not in data:
            continue
        number, err = _as_finite_float(data[key], key, where)
        if err is not None or number is None:
            # 契约：两者恰好其一为 None。两个都判，是为了让后面的比较拿到的
            # ``number`` 确定是 float，而不是把 ``None`` 带进阈值判断。
            errors.append(err or f"{where}: 键 '{key}' 不是可用数值")
            continue
        if number < 0:
            errors.append(f"{where}: 年龄 '{key}' ({number}) 不能为负")
        elif number > MAX_GEOLOGICAL_AGE_MA:
            errors.append(
                f"{where}: 年龄 '{key}' ({number} Ma) 超过地球年龄 "
                f"{MAX_GEOLOGICAL_AGE_MA} Ma"
            )

    for key in spec["prob_keys"]:
        if key not in data:
            continue
        number, err = _as_finite_float(data[key], key, where)
        if err is not None or number is None:
            errors.append(err or f"{where}: 键 '{key}' 不是可用数值")
            continue
        if not 0 < number < 1:
            errors.append(f"{where}: 概率 '{key}' ({number}) 必须落在开区间 (0, 1)")

    for key in spec["positive_keys"]:
        if key not in data:
            continue
        number, err = _as_finite_float(data[key], key, where)
        if err is not None or number is None:
            errors.append(err or f"{where}: 键 '{key}' 不是可用数值")
            continue
        if number <= 0:
            errors.append(f"{where}: '{key}' ({number}) 必须为正值")
        if key == "df" and number <= _MIN_DF:
            errors.append(
                f"{where}: 'df' ({number}) 必须 > {_MIN_DF}，否则 t 分布方差发散"
            )

    # min <= max（区间形状）
    if "min" in spec["required"] and "min" in data and "max" in data:
        lower, _ = _as_finite_float(data["min"], "min", where)
        upper, _ = _as_finite_float(data["max"], "max", where)
        if lower is not None and upper is not None and lower >= upper:
            errors.append(
                f"{where}: 'min' ({lower}) 必须小于 'max' ({upper})；"
                f"上下界写反不会被自动调转（审阅项 §六.33），"
                f"等值区间请用 type: fixed + age"
            )

    note = spec.get("note")
    # 只在用户真的写了那个"一名二义"的键时才附提示，避免对每个正常文件都唠叨
    if note and canonical == "soft_lower" and "p" in data:
        errors.append(f"{where} [{canonical}]: {_NOTE_PREFIX}{note}")

    return errors


def validate_calibration_entry(entry: Any, index: Optional[int] = None) -> List[str]:
    """内容级校验单个校准条目（B-24 的第 1 项修法）。

    Returns:
        信息列表（空表示该条目结构完整、字段合法）。列表里也可能包含
        :data:`_NOTE_PREFIX` 标记的**语义提示**而非错误，调用方用
        :func:`hard_errors` / :func:`informationals` 分流。
    """
    label = f"校准条目 #{index}" if index is not None else "校准条目"
    if not isinstance(entry, dict):
        return [
            f"{label}: 必须是一个映射（- name: ... / constraint: ...），"
            f"实际是 {type(entry).__name__}: {entry!r}"
        ]

    raw_errors: List[str] = []
    if not entry:
        raw_errors.append(f"{label}: 条目为空")

    name = entry.get("name")
    if name is None or (isinstance(name, str) and not name.strip()):
        raw_errors.append(f"{label}: 缺少非空的 'name' 键（收到 {name!r}）")
    elif not isinstance(name, str):
        raw_errors.append(f"{label}: 'name' 必须是字符串，收到 {type(name).__name__}")

    if "constraint" not in entry:
        raw_errors.append(
            f"{label}: 缺少 'constraint' 块（收到键 {sorted(map(str, entry.keys()))}）。"
            f"年龄界必须写在 constraint 之内，不能直接写在条目层"
        )
        constraint_errors: List[str] = []
    else:
        constraint_errors = validate_constraint_block(
            entry["constraint"], where=f"{label} '{name}'.constraint"
        )

    # 节点定义块
    node_label = f"{label} '{name}'.node"
    if "node" in entry:
        raw_errors.extend(_validate_mapping_block(entry["node"], node_label, "node"))
        if isinstance(entry["node"], dict):
            raw_errors.extend(_validate_taxa_and_pair(entry["node"], f"{node_label}"))
    # 顶层 taxa / mrca_pair 与 fossil 块
    raw_errors.extend(_validate_taxa_and_pair(entry, f"{label} '{name}'"))
    if "fossil" in entry:
        fossil_label = f"{label} '{name}'.fossil"
        raw_errors.extend(
            _validate_mapping_block(entry["fossil"], fossil_label, "fossil")
        )
        if isinstance(entry["fossil"], dict):
            raw_errors.extend(_validate_fossil_block(entry["fossil"], fossil_label))
    if "is_root" in entry and not isinstance(entry["is_root"], bool):
        raw_errors.append(
            f"{label} '{name}': 'is_root' 必须是布尔值，收到 "
            f"{type(entry['is_root']).__name__}: {entry['is_root']!r}"
        )

    return raw_errors + constraint_errors


def _validate_fossil_block(block: Dict[str, Any], where: str) -> List[str]:
    """``fossil`` 块的地层年龄：必须是有限数值、非负、min <= max（C-6/B-19）。"""
    errors: List[str] = []
    values: Dict[str, float] = {}
    for key in ("min_age", "max_age"):
        if key not in block or block[key] is None:
            continue
        number, err = _as_finite_float(block[key], key, where)
        if err is not None or number is None:
            errors.append(err or f"{where}: 键 '{key}' 不是可用数值")
            continue
        if number < 0:
            errors.append(f"{where}: '{key}' ({number} Ma) 不能为负")
        elif number > MAX_GEOLOGICAL_AGE_MA:
            errors.append(
                f"{where}: '{key}' ({number} Ma) 超过地球年龄 {MAX_GEOLOGICAL_AGE_MA} Ma"
            )
        else:
            values[key] = number
    if (
        "min_age" in values
        and "max_age" in values
        and values["min_age"] > values["max_age"]
    ):
        errors.append(
            f"{where}: 'min_age' ({values['min_age']}) 不能大于 'max_age' "
            f"({values['max_age']})；两者是含化石地层的'最年轻/最老可能年龄'，"
            f"写反会倒转报告出来的地层区间（不会自动调转）"
        )
    return errors


def _validate_mapping_block(value: Any, where: str, key_name: str) -> List[str]:
    """``node`` / ``fossil`` 这类"必须是映射"的子块。"""
    if value is None:
        return [f"{where}: '{key_name}' 为空（YAML 里写了 '{key_name}:' 但没有内容）"]
    if not isinstance(value, dict):
        return [
            f"{where}: '{key_name}' 必须是一个映射，实际是 "
            f"{type(value).__name__}: {value!r}"
        ]
    return []


def _validate_taxa_and_pair(container: Dict[str, Any], where: str) -> List[str]:
    errors: List[str] = []
    if "taxa" in container:
        taxa = container["taxa"]
        if not isinstance(taxa, list) or not taxa:
            errors.append(
                f"{where}: 'taxa' 必须是非空列表，收到 {type(taxa).__name__}: {taxa!r}"
            )
        elif not all(isinstance(t, str) and t.strip() for t in taxa):
            errors.append(f"{where}: 'taxa' 的每一项都必须是非空字符串: {taxa!r}")
    if "mrca_pair" in container:
        pair = container["mrca_pair"]
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            errors.append(
                f"{where}: 'mrca_pair' 必须是恰好两个叶节点名的列表，收到 {pair!r}"
            )
        elif not all(isinstance(t, str) and t.strip() for t in pair):
            errors.append(f"{where}: 'mrca_pair' 的两项必须是非空字符串: {pair!r}")
    return errors


_NOTE_PREFIX = "提示 — "


def is_informational(message: str) -> bool:
    """校验信息是否为"语义说明"而非真正的错误。"""
    return _NOTE_PREFIX in message


def hard_errors(messages: List[str]) -> List[str]:
    """只保留真正的错误（剔除说明性提示行）。"""
    return [m for m in messages if not is_informational(m)]


def informationals(messages: List[str]) -> List[str]:
    """只保留说明性提示行。"""
    return [m for m in messages if is_informational(m)]


def validate_calibration_payload(data: Any, source: str = "") -> List[str]:
    """内容级校验整份校准文档（YAML 解析结果），返回错误列表。"""
    where = f" ({source})" if source else ""
    if data is None:
        return [f"校准文件为空{where}"]
    if not isinstance(data, dict):
        return [
            f"校准文件顶层必须是映射{where}，实际是 {type(data).__name__}: {data!r}"
        ]
    if "calibrations" not in data:
        return [f"校准文件缺少 'calibrations' 键{where}"]
    entries = data["calibrations"]
    if not isinstance(entries, list):
        return [f"'calibrations' 必须是列表{where}，实际是 {type(entries).__name__}"]
    if not entries:
        return [f"'calibrations' 列表为空{where}：没有任何校准点，定年无法进行"]

    errors: List[str] = []
    for idx, entry in enumerate(entries, 1):
        errors.extend(validate_calibration_entry(entry, idx))
    return errors


class CalibrationLoader:
    """
    校准配置加载器

    从 YAML 文件加载校准点配置
    """

    # 旧/别名约束类型映射（向后兼容）——保留类属性以兼容既有引用，
    # 真身是模块级的 CONSTRAINT_TYPE_ALIASES（校验器与加载器共用同一份表）。
    _CONSTRAINT_TYPE_ALIASES = CONSTRAINT_TYPE_ALIASES

    def __init__(self) -> None:
        self.logger = get_logger()

    def load(self, config_path: Path, strict: bool = False) -> List[CalibrationPoint]:
        """
        从 YAML 文件加载校准点

        Args:
            config_path: YAML 配置文件路径
            strict: 为 True 时，任何被丢弃的条目都会立刻抛
                :class:`~phylodater.core.exceptions.CalibrationError`。
                默认 False：跳过坏条目，但**逐条 error 级留痕 + 汇报
                "请求 N 个 / 生效 M 个"**（审阅项 B-24：丢点不能只留两行 warning）。

        Returns:
            校准点列表
        """
        if not config_path.exists():
            raise FileNotFoundError(f"Calibration config not found: {config_path}")

        with open(config_path, "r") as f:
            try:
                data = yaml.safe_load(f)
            except yaml.YAMLError as e:
                self.logger.error(
                    f"Invalid YAML in calibration config {config_path}: {e}"
                )
                if strict:
                    raise CalibrationError(f"校准配置文件不是合法 YAML: {e}")
                return []

        if not data:
            return []

        entries = data.get("calibrations") or []
        if not isinstance(entries, list):
            self.logger.error(
                f"'calibrations' 必须是列表 {config_path}: 实际是 {type(entries).__name__}"
            )
            if strict:
                raise CalibrationError(f"'calibrations' 必须是列表: {config_path}")
            return []

        requested = len(entries)
        calibrations: List[CalibrationPoint] = []
        dropped: List[str] = []

        # 解析 calibrations 列表（非映射条目也走这里，不再裸 AttributeError，C-37）
        for index, cal_data in enumerate(entries, 1):
            try:
                cal = self._parse_calibration(cal_data, index=index)
            except CalibrationError as e:
                # 内容级校验已经在 _parse_calibration 里 error 过，这里只记账
                dropped.append(f"#{index}: {e}")
                continue
            if cal:
                calibrations.append(cal)
            else:
                dropped.append(f"#{index} 未通过校验（见上面的 error 日志）")

        applied = len(calibrations)
        self.logger.info(
            f"Loaded {applied}/{requested} calibrations from {config_path}"
        )
        if dropped:
            # "N 个点里丢了 k 个" 此前只在日志里留两行 warning（B-24/B-6 同链）
            self.logger.error(
                f"Calibration reconciliation: requested {requested} points but only "
                f"{applied} were applied. {len(dropped)} dropped:\n"
                + "\n".join(f"  - {line}" for line in dropped)
            )
            if strict:
                raise CalibrationError(
                    f"{len(dropped)}/{requested} 个校准条目未能加载，"
                    f"已按 --strict 语义中止: " + "; ".join(dropped)
                )

        return calibrations

    def _parse_calibration(
        self, data: Any, index: Optional[int] = None
    ) -> Optional[CalibrationPoint]:
        """解析单个校准点

        Raises:
            CalibrationError: 条目形状/字段不合法（内容级校验，B-24）。
                加载器**不再**让裸 ``AttributeError`` 冒到用户面前（C-37）。
        """
        label = f"校准条目 #{index}" if index is not None else "校准条目"
        entry_errors = validate_calibration_entry(data, index)
        hard = hard_errors(entry_errors)
        for note in informationals(entry_errors):
            self.logger.warning(note)
        if hard:
            for message in hard:
                self.logger.error(message)
            raise CalibrationError(f"{label} 校验失败（{len(hard)} 处）")

        name = data["name"]

        # 解析年龄约束
        constraint = self._parse_constraint(data.get("constraint"))
        if not constraint:
            self.logger.error(f"Calibration '{name}' has no valid constraint, skipping")
            return None

        # 解析节点定义（支持顶层 taxa/mrca_pair 或 node 块内定义）
        resolved_taxa = data.get("taxa") or []
        mrca_leaf_pair = tuple(data["mrca_pair"]) if data.get("mrca_pair") else None

        node_def = None
        if "node" in data:
            try:
                node_def = NodeDefinition.from_dict(data["node"])
            except (TypeError, ValueError, KeyError, AttributeError) as e:
                raise CalibrationError(
                    f"{label} '{name}'.node: 无法解析节点定义 ({e})"
                ) from e
            # node 块内定义的 taxa / mrca_pair 优先级更高
            if node_def.taxa:
                resolved_taxa = node_def.taxa
            if node_def.mrca_pair:
                mrca_leaf_pair = tuple(node_def.mrca_pair)

        # 解析化石元数据
        fossil_metadata = None
        if "fossil" in data:
            try:
                fossil_metadata = self._parse_fossil_metadata(data["fossil"])
            except (TypeError, ValueError) as e:
                raise CalibrationError(
                    f"{label} '{name}'.fossil: 化石元数据不合法 ({e})"
                ) from e

        try:
            return CalibrationPoint(
                name=name,
                age_constraint=constraint,
                resolved_taxa=resolved_taxa,
                mrca_leaf_pair=mrca_leaf_pair,
                is_root_node=bool(data.get("is_root", False)),
                fossil_metadata=fossil_metadata,
            )
        except (TypeError, ValueError) as e:
            raise CalibrationError(
                f"{label} '{name}': 无法构造 CalibrationPoint ({e})"
            ) from e

    def _parse_constraint(self, data: Optional[Dict[str, Any]]) -> Optional[Any]:
        """解析年龄约束（内容级校验见 :func:`validate_constraint_block`）

        审阅项 B-19 / B-24 / C-38 带来的三点改变：

        1. **不再猜类型**：缺 ``type`` 直接判为无效并说明实际收到的键，而不是
           静默按 ``fixed`` 处理、再报一句指错方向的 "缺少 age"；
        2. **非有限年龄在建构前就被拒**：YAML 的 ``.nan`` / ``.inf`` 字面量会被
           pyyaml 解析成 ``float('nan')`` / ``float('inf')``，而
           ``float()`` 之前原来什么检查都没有；
        3. **失败一律 error 级留痕**：``except (KeyError, ValueError)`` 把问题压成
           一行 warning 然后丢点，是 B-24 里"校准静默缩水"的直接成因。
        """
        messages = validate_constraint_block(data)
        hard = hard_errors(messages)
        for note in informationals(messages):
            self.logger.warning(note)
        # ``validate_constraint_block`` 对非映射、缺 type、未知 type 一律报硬错误，
        # 所以这里的早退已经保证 ``data`` 是可下标的字典。把收窄结果存进
        # ``block``：局部变量 ``data`` 的收窄不会传递到嵌套函数 ``num`` 里。
        if hard or not isinstance(data, dict):
            for message in hard:
                self.logger.error(message)
            return None
        block: Dict[str, Any] = data

        constraint_type = _canonical_type_name(block["type"])
        if constraint_type in CONSTRAINT_TYPE_ALIASES:
            self.logger.warning(
                f"Deprecated constraint type '{block['type']}' translated to "
                f"'{constraint_type}'"
            )

        def num(key: str, default: Optional[float] = None) -> float:
            if key in block and block[key] is not None:
                value, err = _as_finite_float(block[key], key, "constraint")
                if err is not None or value is None:  # pragma: no cover - 上层已拦
                    raise ValueError(err or f"constraint: 键 '{key}' 不是可用数值")
                return value
            if default is None:
                raise KeyError(key)
            return default

        try:
            if constraint_type == "fixed":
                return FixedAgeConstraint(fixed_age=num("age"))

            elif constraint_type == "uniform":
                return UniformAgeConstraint(
                    min_age=num("min"),
                    max_age=num("max"),
                    **(
                        {
                            "tail_lower": num("p_lower"),
                            "tail_upper": num("p_upper"),
                        }
                        if ("p_lower" in block or "p_upper" in block)
                        else {}
                    ),
                )

            elif constraint_type == "soft_lower":
                # YAML 的 p/c 就是 PAML 手册 L(tL, p, c) 的 p 与 c：
                # p = 相对偏移比例 offset_fraction，c = Cauchy 尺度 cauchy_scale。
                # （旧代码写的是 tail_prob=/tail_shape=，字段已随审阅项 B-2 改名。）
                kwargs: Dict[str, float] = {
                    "min_age": num("min"),
                    "offset_fraction": num("p", 0.1),
                    "cauchy_scale": num("c", 1.0),
                }
                if "tail_prob" in block:
                    # L() 的第 4 槽（真正的左尾概率），旧版本从不填
                    kwargs["tail_prob"] = num("tail_prob")
                return SoftLowerBoundConstraint(**kwargs)

            elif constraint_type == "maximum":
                return MaximumAgeConstraint(
                    max_age=num("max"), tail_prob=num("p", 0.025)
                )

            elif constraint_type == "soft_bounds":
                if "p" in block and ("p_lower" in block or "p_upper" in block):
                    raise ValueError(
                        "soft_bounds: 'p'（对称尾部概率）与 'p_lower'/'p_upper' 不能同时给出"
                    )
                soft_kwargs: Dict[str, float] = {
                    "min_age": num("min"),
                    "max_age": num("max"),
                }
                if "p" in block:
                    soft_kwargs["tail_prob"] = num("p")
                if "p_lower" in block:
                    soft_kwargs["tail_lower"] = num("p_lower")
                if "p_upper" in block:
                    soft_kwargs["tail_upper"] = num("p_upper")
                return SoftBoundsConstraint(**soft_kwargs)

            elif constraint_type == "gamma":
                return GammaPriorConstraint(
                    alpha=num("alpha"),
                    beta=num("beta"),
                    offset=num("offset", 0.0),
                    scale=num("scale", 1.0),
                )

            elif constraint_type == "skew_normal":
                return SkewNormalConstraint(
                    location=num("location"),
                    scale=num("scale"),
                    shape=num("shape"),
                )

            elif constraint_type == "skew_t":
                return SkewTConstraint(
                    location=num("location"),
                    scale=num("scale"),
                    shape=num("shape"),
                    df=num("df"),
                )

            else:  # pragma: no cover - validate_constraint_block 已拦
                self.logger.error(f"Unknown constraint type: {block.get('type')!r}")
                return None

        except (KeyError, TypeError, ValueError) as e:
            # 约束构造期校验（constraints.py）也可能拒绝：那必须让用户看见
            self.logger.error(
                f"Failed to build {constraint_type} constraint: {e} "
                f"(received keys: {sorted(map(str, block.keys()))})"
            )
            return None

    def _parse_fossil_metadata(self, data: Dict[str, Any]) -> FossilMetadata:
        """解析化石元数据（非有限年龄由 :class:`FossilMetadata` 拒绝）"""
        confidence = None
        if "confidence" in data:
            try:
                confidence = FossilConfidence(data["confidence"])
            except ValueError:
                self.logger.warning(
                    f"化石元数据的 confidence '{data['confidence']}' 不是 "
                    f"confirmed/probable/tentative 之一，已忽略该字段"
                )

        return FossilMetadata(
            fossil_name=data.get("name"),
            formation=data.get("formation"),
            geological_period=data.get("period"),
            min_stratigraphic_age=data.get("min_age"),
            max_stratigraphic_age=data.get("max_age"),
            reference=data.get("reference"),
            confidence=confidence,
            notes=data.get("notes"),
        )

    def save_example_config(self, output_path: Path) -> None:
        """生成示例配置文件"""
        example = {
            "calibrations": [
                {
                    "name": "LUCA",
                    "constraint": {"type": "uniform", "min": 3500, "max": 4500},
                    "node": {"type": "auto"},
                    "is_root": True,
                    "fossil": {
                        "name": "Earliest life evidence",
                        "period": "Archean",
                        "reference": "doi:10.xxxx/xxxxx",
                    },
                },
                {
                    "name": "Cyanobacteriota",
                    "constraint": {
                        "type": "soft_lower",
                        "min": 2000,
                        "p": 0.1,
                        "c": 1.0,
                    },
                    "node": {"type": "auto", "taxa": ["Cyanobacteriota"]},
                },
                {
                    "name": "Proteobacteria",
                    "constraint": {"type": "fixed", "age": 1500},
                    "node": {
                        "type": "mrca",
                        "taxa": ["Alphaproteobacteria", "Gammaproteobacteria"],
                    },
                },
            ]
        }

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with safe_writer(output_path) as f:
            yaml.dump(example, f, default_flow_style=False, sort_keys=False)

        self.logger.info(f"Generated example config: {output_path}")


def create_default_calibration_config(output_path: Path) -> None:
    """创建默认校准配置文件"""
    loader = CalibrationLoader()
    loader.save_example_config(output_path)
