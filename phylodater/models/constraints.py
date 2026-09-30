"""
年龄约束类型体系

定义所有年龄约束类型及其向各软件格式的转换
"""

import math
import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from scipy import stats

MAX_GEOLOGICAL_AGE_MA = 4600.0
MAX_GEOLOGICAL_AGE_GA = MAX_GEOLOGICAL_AGE_MA / 1000.0

# --- MCMCTree 侧的先验参数默认值（与上游逐字对应） ------------------------------- #
# treesub.c:8675:
#     double tailL = 0.025, tailR = 0.025, p_LOWERBOUND = 0.1, c_LOWERBOUND = 1.0;
# 这些值只在注释串"没有填满槽位"时生效；本项目一律填满槽位，因此它们只作为
# "与上游一致的默认值"被显式写出（审阅项 B-1/B-2）。
MCMCTREE_DEFAULT_TAIL_PROB = 0.025
MCMCTREE_DEFAULT_LOWER_P = 0.1
MCMCTREE_DEFAULT_LOWER_C = 1.0

# 固定年龄用 ±1% 的窄软边界近似（MCMCTree 无点校准）；每侧 0.5% 尾部概率
# ⇒ 99% 的先验概率质量落在 fixed_age ±1% 窗口内（审阅项 A-3/B-1）。
FIXED_AGE_TOLERANCE_FRACTION = 0.01
FIXED_AGE_TAIL_PROB = 0.005

# MCMCTree 也**没有**真正的硬均匀界：BOUND_F 密度把 tailL/tailR 用作除数
# （mcmctree.c:2613-2620，thetaL = (1-tailL-tailR)*a/(tailL*(b-a))），传 0 会在
# t<a / t>b 分支上得到 0*inf = NaN。因此"硬界"只能用极小的非零尾部近似。
NEAR_HARD_TAIL_PROB = 0.0001

# MCMCTree 串的打印精度（见 `_format_ga`：默认 4 位小数，不足时自动扩展）
_MCMCTREE_GA_DECIMALS = 4
_MCMCTREE_GA_SIG_DIGITS = 6
_MCMCTREE_NUMBER_DECIMALS = 10


def _emit_degradation_warning(msg: str) -> None:
    """发出语义降级警告（非致命）"""
    from ..core.exceptions import SemanticDegradationWarning

    warnings.warn(msg, SemanticDegradationWarning, stacklevel=3)


def _require_finite_number(value: object, name: str) -> float:
    """类型 + 有限性校验（审阅项 B-19）。

    MCMCTree 侧所有越界检查都用 ``>`` / ``<=`` 比较，而 NaN 与任何数比较均为
    False，所以 NaN 年龄能穿过全部数值校验、并被 ``sscanf("%lf")`` 原样读入，
    最终让先验对数密度在整个定义域上恒为 NaN。inf 过去会被
    ``MAX_GEOLOGICAL_AGE_MA`` 守卫拦下，因此这里把两者一起在建构期拒掉。
    """
    if not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a number, got {type(value).__name__}")
    if not math.isfinite(value):
        raise ValueError(
            f"{name} ({value}) must be a finite real number: NaN and infinity are "
            f"rejected because every bound check uses >/< comparisons that NaN "
            f"silently passes"
        )
    return float(value)


def _format_number(value: float) -> str:
    """把概率/形状参数渲染为定点小数字符串（不产出 ``1e-05`` 这类科学计数）。"""
    text = f"{value:.{_MCMCTREE_NUMBER_DECIMALS}f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def _format_ga(ga: float) -> str:
    """Ga 数值 → 字符串。

    默认沿用 4 位小数（与历史输出/文档一致）；当 4 位小数会把两个不同的界打印成
    同一个数（即 ``float()`` 回读不相等）时，按有效数字自动扩展小数位。
    这是审阅项 A-3 的一半修复：0.1/0.5/1.0 Ma 级别的校准不再坍缩为零宽度。
    """
    text = f"{ga:.{_MCMCTREE_GA_DECIMALS}f}"
    if float(text) == ga:
        return text
    exponent = math.floor(math.log10(abs(ga)))
    decimals = max(_MCMCTREE_GA_DECIMALS, _MCMCTREE_GA_SIG_DIGITS - 1 - exponent)
    text = f"{ga:.{decimals}f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def _render_ga_bound(age_ma: float, name: str, source: str) -> str:
    """Ma → Ga 字符串，并校验**最终渲染值**落在 [0, 地球年龄]（审阅项 C-2/B-19）。

    输入侧的 ``MAX_GEOLOGICAL_AGE_MA`` 守卫只看用户写进的年龄；软化和 ±1% 容差
    窗口会把上界推过 4.6 Ga（例如 ``fixed_age=4600`` → 4.646 Ga），因此必须在
    出口处再校验一次。
    """
    if not math.isfinite(age_ma):
        raise ValueError(f"{source}: {name} ({age_ma} Ma) is not a finite number")
    ga = age_ma / 1000.0
    if ga < 0:
        raise ValueError(f"{source}: {name} ({age_ma} Ma) must be non-negative")
    if ga > MAX_GEOLOGICAL_AGE_GA * (1.0 + 1e-9):
        raise ValueError(
            f"{source}: rendered {name} = {_format_ga(ga)} Ga ({age_ma} Ma) exceeds "
            f"the maximum geological age ({MAX_GEOLOGICAL_AGE_MA} Ma = "
            f"{MAX_GEOLOGICAL_AGE_GA} Ga, Earth's age). Narrow the calibration "
            f"window (a ±1% soft window requires fixed_age <="
            f" {MAX_GEOLOGICAL_AGE_MA / (1 + FIXED_AGE_TOLERANCE_FRACTION):.2f} Ma)."
        )
    return _format_ga(ga)


def _validate_tail_probs(source: str, tail_lower: float, tail_upper: float) -> None:
    """校验 BOUND_F 的两个尾部概率槽（审阅项 B-1/B-2）。"""
    tail_lower = _require_finite_number(tail_lower, f"{source}.tail_lower")
    tail_upper = _require_finite_number(tail_upper, f"{source}.tail_upper")
    for name, value in (("tail_lower", tail_lower), ("tail_upper", tail_upper)):
        if not 0 < value < 1:
            raise ValueError(
                f"{source}: {name} ({value}) must be a probability in (0, 1); "
                f"MCMCTree uses it as a divisor in the BOUND_F density, so 0 is not "
                f"a legal 'hard bound' switch"
            )
    if tail_lower + tail_upper >= 1:
        raise ValueError(
            f"{source}: tail_lower + tail_upper ({tail_lower + tail_upper}) must be "
            f"< 1, otherwise the flat part of the BOUND_F prior carries no probability"
        )


def _render_mcmctree_bounds(
    source: str,
    lower_ma: float,
    upper_ma: float,
    tail_lower: float,
    tail_upper: float,
) -> str:
    """渲染 MCMCTree 的 ``B(tL, tU, tailL, tailR)``（4 个槽，审阅项 B-1）。

    上游期望 4 个值（``npfossils[BOUND_F] = 4``，treesub.c:8744-8748）：只写 3 个
    时第 3 槽被当作 **tailL**，tailR 保持栈上默认 0.025，于是"对称软界"实际是非
    对称的，而"硬界"（p=0）会把 tailL 置 0 并在 t<a 分支上除零。
    """
    lower_text = _render_ga_bound(lower_ma, "lower bound", source)
    upper_text = _render_ga_bound(upper_ma, "upper bound", source)
    _validate_tail_probs(source, tail_lower, tail_upper)
    lower_ga, upper_ga = float(lower_text), float(upper_text)
    if not upper_ga > lower_ga:
        raise ValueError(
            f"{source}: rendered MCMCTree interval B({lower_text}, {upper_text}) is "
            f"degenerate (zero width) or inverted (lower >= upper). MCMCTree aborts on "
            f"lower > upper ('fossil bounds in tree incorrect', treesub.c:8749-8752) "
            f"and divides by (b - a) in the BOUND_F density (mcmctree.c:2611-2619), "
            f"so this prior is unusable as rendered. Coarsen the formatting or widen "
            f"the interval."
        )
    return (
        f"B({lower_text}, {upper_text}, "
        f"{_format_number(tail_lower)}, {_format_number(tail_upper)})"
    )


def _render_mcmctree_upper_bound(
    source: str, upper_ma: float, tail_upper: float
) -> str:
    """渲染 ``U(tU, tailR)``：上游只读 2 个槽（treesub.c:8740-8743）。"""
    upper_text = _render_ga_bound(upper_ma, "upper bound", source)
    tail_upper = _require_finite_number(tail_upper, f"{source}.tail_prob")
    if not 0 < tail_upper < 1:
        raise ValueError(
            f"{source}: tail_prob ({tail_upper}) must be a probability in (0, 1)"
        )
    return f"U({upper_text}, {_format_number(tail_upper)})"


def _render_mcmctree_lower_bound(
    source: str,
    lower_ma: float,
    offset_fraction: float,
    cauchy_scale: float,
    tail_lower: float,
) -> str:
    """渲染 ``L(tL, P, c, tailL)``（4 个槽，审阅项 B-2）。

    上游槽位语义（treesub.c:8732-8738 + mcmctree.c:2581-2596）：
    ``P`` 是先验中心的**相对偏移比例**（``t0 = tL*(1+P)``）、``c`` 是截断 Cauchy
    的尺度（``s = tL*c``）、第 4 槽才是左尾概率 ``tailL``。只写 3 个值时 tailL
    永远停在默认 0.025，用户设的尾概率到不了那里。
    """
    lower_text = _render_ga_bound(lower_ma, "lower bound", source)
    offset_fraction = _require_finite_number(
        offset_fraction, f"{source}.offset_fraction"
    )
    cauchy_scale = _require_finite_number(cauchy_scale, f"{source}.cauchy_scale")
    tail_lower = _require_finite_number(tail_lower, f"{source}.tail_prob")
    if offset_fraction < 0:
        raise ValueError(
            f"{source}: offset_fraction ({offset_fraction}) must be >= 0; it shifts "
            f"the Cauchy centre to tL*(1+offset_fraction)"
        )
    if cauchy_scale <= 0:
        raise ValueError(f"{source}: cauchy_scale ({cauchy_scale}) must be positive")
    if not 0 < tail_lower < 1:
        raise ValueError(
            f"{source}: tail_prob ({tail_lower}) must be a probability in (0, 1)"
        )
    return (
        f"L({lower_text}, {_format_number(offset_fraction)}, "
        f"{_format_number(cauchy_scale)}, {_format_number(tail_lower)})"
    )


def _validate_derived_interval(
    source: str, ci_lower: float, ci_upper: float
) -> Tuple[float, float, List[str]]:
    """校验由先验矩/分位数推导出的降级区间（审阅项 B-19/C-2/C-5）。

    返回 ``(ci_lower, ci_upper, notes)``。``notes`` 记录为让区间可写而做的**语义改动**
    （目前是把负下界截断到 0——这使区间不再是等尾 95% 区间）；调用方必须用
    :func:`_notes_suffix` 把它写进降级警告文本，否则警告会描述一个并不是真正写出去
    的区间（C-5 的原始缺陷正是"静默截断、只字不提"）。
    """
    notes: List[str] = []
    ci_lower = _require_finite_number(ci_lower, f"{source}: derived lower bound")
    ci_upper = _require_finite_number(ci_upper, f"{source}: derived upper bound")
    if not ci_upper > ci_lower:
        raise ValueError(
            f"{source}: derived degradation interval [{ci_lower}, {ci_upper}] Ma is "
            f"degenerate (zero width) or inverted; refusing to render it"
        )
    if ci_lower < 0:
        notes.append(
            f"the lower bound was truncated from {ci_lower:.2f} Ma to 0.0 Ma, so this "
            f"is no longer an equal-tailed 95% interval (probability mass below 0 is "
            f"cut off rather than folded back)"
        )
        ci_lower = 0.0
    if ci_upper > MAX_GEOLOGICAL_AGE_MA:
        notes.append(
            f"the upper bound {ci_upper:.1f} Ma exceeds the maximum geological age "
            f"({MAX_GEOLOGICAL_AGE_MA} Ma)"
        )
    return ci_lower, ci_upper, notes


def _notes_suffix(notes: Optional[Sequence[str]]) -> str:
    """把 :func:`_validate_derived_interval` 的改动清单拼成警告后缀。"""
    if not notes:
        return ""
    return " NOTE: " + " ".join(notes)


class AgeConstraint(ABC):
    """年龄约束抽象基类"""

    @abstractmethod
    def _render_software_format(self, target_software: str) -> str:
        """构建目标软件格式的约束模板（可含 {label}/{Name}/{t1}/{t2}/{taxa_str} 占位符）。

        由 :meth:`to_software_format` 代入节点真实字段后返回最终字符串。
        """
        pass

    def to_software_format(
        self,
        target_software: str,
        *,
        label: str = "",
        name: str = "",
        t1: str = "",
        t2: str = "",
        taxa_str: str = "",
    ) -> str:
        """返回目标软件格式、已代入节点真实字段的约束字符串（不含占位符）。

        各适配器在调用时传入与模板占位符对应的节点字段，模型内部完成渲染，
        避免适配器自行 ``.format(...)`` 时因占位符键名随约束分支不同而 KeyError
        （旧实现的脆弱点）。

        Args:
            target_software: 目标软件名（mcmctree/r8s/pyr8s/treepl/pathd8/lsd2/
                wlogdate/mdcat）。
            label: r8s/pyr8s 使用的节点标签（对应模板 ``{label}``）。
            name: treepl/wlogdate/mdcat 使用的校准点名（对应模板 ``{Name}``）。
            t1, t2: pathd8 使用的 MRCA 两叶名（对应模板 ``{t1}``/``{t2}``）。
            taxa_str: wlogdate/mdcat 使用的叶节点串（'+' 连接，对应模板
                ``{taxa_str}``）。
        """
        template = self._render_software_format(target_software)
        return template.format(label=label, Name=name, t1=t1, t2=t2, taxa_str=taxa_str)

    @abstractmethod
    def to_mcmctree_calib_string(self) -> str:
        """转换为 MCMCTree 校准字符串"""
        pass

    def describe(self) -> str:
        """返回约束的人类可读摘要（用于报告/日志/校准树文件）。

        子类应覆盖此方法，提供完整的约束描述；默认实现仅返回类名。
        """
        return f"{type(self).__name__}"


@dataclass
class FixedAgeConstraint(AgeConstraint):
    """固定年龄约束（MCMCTree 侧以 ±1% 窄软界 ``B()`` 近似，见
    :meth:`to_mcmctree_calib_string`）"""

    fixed_age: float  # 单位：Ma

    def __post_init__(self) -> None:
        _require_finite_number(self.fixed_age, "fixed_age")
        if self.fixed_age <= 0:
            raise ValueError(f"fixed_age ({self.fixed_age}) must be positive")
        if self.fixed_age > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"fixed_age ({self.fixed_age} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma - Earth's age)"
            )

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        elif target_software in ["r8s", "pyr8s"]:
            return f"fixage taxon={{label}} age={self.fixed_age};"
        elif target_software == "treepl":
            # treePL 使用等值的 min 和 max
            return f"min = {{Name}} {self.fixed_age}\nmax = {{Name}} {self.fixed_age}"
        elif target_software == "pathd8":
            return f"mrca: {{t1}}, {{t2}}, fixage={self.fixed_age};"
        elif target_software == "lsd2":
            # LSD2 使用负值
            return f"-{self.fixed_age}"
        elif target_software == "wlogdate":
            return f"{{Name}}={{taxa_str}}\t{self.fixed_age}"
        elif target_software == "mdcat":
            return f"{{Name}}={{taxa_str}}\t{self.fixed_age}"
        else:
            raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """把固定年龄近似成 MCMCTree 的窄软边界 ``B(tL, tU, tailL, tailR)``。

        MCMCTree 不支持点校准（``@X`` 语法只属于 r8s/PATHd8 一类软件），只能给
        BOUND_F 先验；本实现用 ``fixed_age ± 1%`` 的窗口、每侧 0.5% 尾部概率，
        使 99% 的先验概率质量落在 ±1% 窗口内。

        审阅项 A-3/C-2/B-1 的三点修正：

        1. **精度与容差解耦**：旧实现固定用 ``.4f``（可分辨步长 1e-4 Ga = 0.1 Ma），
           而容差只有年龄的 1%（= ``age_Ma * 1e-5`` Ga），于是 age < 5 Ma 时上下界
           被打印成同一个数（0.1/0.5/1.0 Ma → 零宽度）。现由 :func:`_format_ga`
           按需扩展小数位，渲染后严格 ``lower < upper``。
        2. **不再单向夹紧**：旧实现 ``min_age = max(0.0001, age_ga - tolerance)`` 只抬
           下界、不抬上界，age ≲ 0.1 Ma 时得到 ``B(0.0001, 0.0000)``（下界 > 上界）。
           由于容差恒为 1%，``0.99 * age`` 永远为正、永远低于上界，夹紧地板本身不再
           需要；万一渲染结果退化则由 :func:`_render_mcmctree_bounds` 抛错。
        3. **4 个槽写满**：旧串只给 3 个值，第 3 值被上游当作 tailL、tailR 停在默认
           0.025（不对称软界）。现显式给出 tailL/tailR。
        """
        tolerance = self.fixed_age * FIXED_AGE_TOLERANCE_FRACTION
        return _render_mcmctree_bounds(
            source=f"FixedAgeConstraint(fixed_age={self.fixed_age} Ma)",
            lower_ma=self.fixed_age - tolerance,
            upper_ma=self.fixed_age + tolerance,
            tail_lower=FIXED_AGE_TAIL_PROB,
            tail_upper=FIXED_AGE_TAIL_PROB,
        )

    def describe(self) -> str:
        """人类可读摘要：固定年龄点估计（MCMCTree 侧以 ±1% 窄软界近似）。"""
        return (
            f"{self.fixed_age} Ma (fixed point; MCMCTree approximates it as "
            f"±{FIXED_AGE_TOLERANCE_FRACTION * 100:.0f}% soft bounds)"
        )


@dataclass
class UniformAgeConstraint(AgeConstraint):
    """均匀分布年龄约束（硬边界，MCMCTree 侧以极小尾部的 BOUND_F 近似）"""

    min_age: float  # 单位：Ma
    max_age: float  # 单位：Ma
    # MCMCTree BOUND_F 的尾部概率槽（第 3/4 槽）。上游没有"硬界"开关：0 会作为除数
    # 出现在 mcmctree.c:2613-2620，因此默认取极小但非零的值（每侧 0.01%，即 99.98%
    # 概率质量在区间内），需要更松的界时把它调到 0.025（上游默认值）。
    tail_lower: float = NEAR_HARD_TAIL_PROB
    tail_upper: float = NEAR_HARD_TAIL_PROB

    def __post_init__(self) -> None:
        _require_finite_number(self.min_age, "min_age")
        _require_finite_number(self.max_age, "max_age")
        if self.min_age <= 0:
            raise ValueError(f"min_age ({self.min_age}) must be positive")
        if self.max_age > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"max_age ({self.max_age} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma - Earth's age)"
            )
        if self.min_age >= self.max_age:
            raise ValueError(
                f"min_age ({self.min_age}) must be < max_age ({self.max_age})"
            )

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        elif target_software in ["r8s", "pyr8s"]:
            return (
                f"constrain taxon={{label}} min_age={self.min_age};\n"
                f"constrain taxon={{label}} max_age={self.max_age};"
            )
        elif target_software == "treepl":
            return f"min = {{Name}} {self.min_age}\n" f"max = {{Name}} {self.max_age}"
        elif target_software == "pathd8":
            # PATHd8 必须分写为两行
            return (
                f"mrca: {{t1}}, {{t2}}, minage={self.min_age};\n"
                f"mrca: {{t1}}, {{t2}}, maxage={self.max_age};"
            )
        elif target_software == "lsd2":
            # LSD2 使用负值，绝对值较大者在前
            return f"-{self.max_age}:-{self.min_age}"
        elif target_software == "wlogdate":
            # wLogDate 仅支持点估计，降级为区间中值
            mid = (self.min_age + self.max_age) / 2
            _emit_degradation_warning(
                f"wLogDate does not support range constraints. "
                f"Using midpoint {mid} from range [{self.min_age}, {self.max_age}]"
            )
            return f"{{Name}}={{taxa_str}}\t{mid}"
        elif target_software == "mdcat":
            # MD-Cat 仅支持点估计，降级为区间中值
            mid = (self.min_age + self.max_age) / 2
            _emit_degradation_warning(
                f"MD-Cat does not support range constraints. "
                f"Using midpoint {mid} from range [{self.min_age}, {self.max_age}]"
            )
            return f"{{Name}}={{taxa_str}}\t{mid}"
        else:
            raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """MCMCTree ``B(tL, tU, tailL, tailR)``（单位 Ga，4 个槽写满）。

        审阅项 B-1：旧串是 ``B(tL, tU, 0.0)``，注释声称"p=0 表示硬边界"，但上游把
        第 3 槽当 **tailL**（mcmctree.c:2609-2610），传 0 会让 ``t < tL`` 分支
        ``thetaL = (1-tailL-tailR)*a/(tailL*(b-a))`` 除零、``log(tailL*thetaL/a)``
        成为 ``log(0·∞)=NaN``，而 tailR 仍悄悄停在默认 0.025。现在两侧尾部都显式
        给出，且取极小非零值来逼近"硬界"（MCMCTree 无真正的硬均匀界）。
        """
        return _render_mcmctree_bounds(
            source=(
                f"UniformAgeConstraint(min_age={self.min_age} Ma, "
                f"max_age={self.max_age} Ma)"
            ),
            lower_ma=self.min_age,
            upper_ma=self.max_age,
            tail_lower=self.tail_lower,
            tail_upper=self.tail_upper,
        )

    def describe(self) -> str:
        """人类可读摘要：均匀区间约束（近似硬界）。"""
        return (
            f"{self.min_age}-{self.max_age} Ma (uniform, near-hard bounds; "
            f"tails={self.tail_lower}/{self.tail_upper})"
        )


@dataclass(init=False)
class SoftLowerBoundConstraint(AgeConstraint):
    """软下界约束（MCMCTree 截断 Cauchy 先验 ``L(tL, P, c, tailL)``）

    字段名与上游槽位一一对应（审阅项 B-2）：

    * ``offset_fraction`` → 槽 2 ``P``：先验中心的**相对偏移比例**，
      ``t0 = tL * (1 + P)``（mcmctree.c:2585）。旧字段名 ``tail_prob`` 把这个量
      错描述成"尾概率"，导致用户以为在收紧尾部、实际在移动先验中心。
    * ``cauchy_scale`` → 槽 3 ``c``：Cauchy 尺度，``s = tL * c``（旧名 ``tail_shape``）。
    * ``tail_prob`` → 槽 4 ``tailL``：**真正的**左尾概率（t < tL 的概率质量）。
      旧实现从不填这一槽，它永远停在上游默认 0.025。
    """

    min_age: float  # 最小年龄（Ma），对应上游 tL
    offset_fraction: float = MCMCTREE_DEFAULT_LOWER_P  # 槽 2：P（相对偏移比例）
    cauchy_scale: float = MCMCTREE_DEFAULT_LOWER_C  # 槽 3：c（Cauchy 尺度）
    tail_prob: float = MCMCTREE_DEFAULT_TAIL_PROB  # 槽 4：tailL（左尾概率）

    def __init__(
        self,
        min_age: float,
        offset_fraction: Optional[float] = None,
        cauchy_scale: Optional[float] = None,
        tail_prob: Optional[float] = None,
        tail_shape: Optional[float] = None,
    ) -> None:
        """兼容旧字段名的一次性换算（审阅项 B-2 修法：改名但不静默改语义）。

        旧 API 为 ``(min_age, tail_prob, tail_shape)``，其中 ``tail_prob`` 实际是
        上游的 ``P``。只有 ``tail_shape`` 这个只存在于旧 API 的关键字能无歧义地
        识别旧调用：此时把 ``tail_prob`` 读作 ``offset_fraction``、``tail_shape``
        读作 ``cauchy_scale``，并发出 DeprecationWarning。
        """
        if tail_shape is not None:
            if offset_fraction is not None:
                raise TypeError(
                    "SoftLowerBoundConstraint: cannot combine legacy 'tail_shape' "
                    "with new 'offset_fraction'; use offset_fraction/cauchy_scale/"
                    "tail_prob only"
                )
            warnings.warn(
                "SoftLowerBoundConstraint(tail_prob=..., tail_shape=...) is "
                "deprecated: 'tail_prob' used to be placed in MCMCTree's L() slot 2, "
                "which is the relative offset P (t0 = tL*(1+P)), not a probability. "
                "Renamed to offset_fraction/cauchy_scale; the numeric value you pass "
                "keeps the same meaning (P -> offset_fraction, c -> tail_shape). "
                "The new tail_prob field maps to slot 4 (tailL).",
                DeprecationWarning,
                stacklevel=2,
            )
            offset_fraction = tail_prob if tail_prob is not None else None
            cauchy_scale = tail_shape
            tail_prob = None

        self.min_age = min_age
        self.offset_fraction = (
            MCMCTREE_DEFAULT_LOWER_P if offset_fraction is None else offset_fraction
        )
        self.cauchy_scale = (
            MCMCTREE_DEFAULT_LOWER_C if cauchy_scale is None else cauchy_scale
        )
        self.tail_prob = MCMCTREE_DEFAULT_TAIL_PROB if tail_prob is None else tail_prob
        self.__post_init__()

    def __post_init__(self) -> None:
        _require_finite_number(self.min_age, "min_age")
        _require_finite_number(self.offset_fraction, "offset_fraction")
        _require_finite_number(self.cauchy_scale, "cauchy_scale")
        _require_finite_number(self.tail_prob, "tail_prob")
        if self.min_age <= 0:
            raise ValueError(f"min_age ({self.min_age}) must be positive")
        if self.min_age > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"min_age ({self.min_age} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma)"
            )
        if self.offset_fraction < 0:
            raise ValueError(
                f"offset_fraction ({self.offset_fraction}) must be non-negative"
            )
        if self.cauchy_scale <= 0:
            raise ValueError(f"cauchy_scale ({self.cauchy_scale}) must be positive")
        if not 0 < self.tail_prob < 1:
            raise ValueError(
                f"tail_prob ({self.tail_prob}) must be a probability strictly "
                f"between 0 and 1 (MCMCTree slot 4 of L())"
            )

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        elif target_software == "lsd2":
            # LSD2 date 文件取值轴：value = -age（负值表示距今），``v1:v2`` 表示
            # value ∈ [v1, v2]（NA 表示该侧无界）。经 IQ-TREE 2.3.6 / 3.1.3 实测
            # 回读校准节点年龄确认：``NA:-T`` 被 LSD 解码为“年龄 ≥ T”（下界），
            # ``-T:NA`` 解码为“年龄 ≤ T”（上界）。
            # 下界约束：年龄至少 min_age → 上端（更老侧）开放 -> NA:-min_age
            return f"NA:-{self.min_age}"
        else:
            # 其他软件降级为硬下界
            _emit_degradation_warning(
                f"{target_software} does not support soft lower bounds. "
                f"Degrading to hard minimum age {self.min_age}"
            )
            if target_software in ["r8s", "pyr8s"]:
                return f"constrain taxon={{label}} min_age={self.min_age};"
            elif target_software == "treepl":
                return f"min = {{Name}} {self.min_age}"
            elif target_software == "pathd8":
                return f"mrca: {{t1}}, {{t2}}, minage={self.min_age};"
            elif target_software in ["wlogdate", "mdcat"]:
                return f"{{Name}}={{taxa_str}}\t{self.min_age}"
            else:
                raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """MCMCTree ``L(tL, P, c, tailL)``（单位 Ga，4 个槽写满，审阅项 B-2）。"""
        return _render_mcmctree_lower_bound(
            source=f"SoftLowerBoundConstraint(min_age={self.min_age} Ma)",
            lower_ma=self.min_age,
            offset_fraction=self.offset_fraction,
            cauchy_scale=self.cauchy_scale,
            tail_lower=self.tail_prob,
        )

    def describe(self) -> str:
        """人类可读摘要：软下界约束（截断 Cauchy）。"""
        return (
            f"≥{self.min_age} Ma (soft lower bound, truncated Cauchy: "
            f"P={self.offset_fraction}, c={self.cauchy_scale}, "
            f"tailL={self.tail_prob})"
        )


@dataclass
class MaximumAgeConstraint(AgeConstraint):
    """最大年龄约束（上界）"""

    max_age: float  # 最大年龄（Ma）
    tail_prob: float = MCMCTREE_DEFAULT_TAIL_PROB  # 尾部概率（上游 tailR，默认 0.025）

    def __post_init__(self) -> None:
        _require_finite_number(self.max_age, "max_age")
        _require_finite_number(self.tail_prob, "tail_prob")
        if self.max_age <= 0:
            raise ValueError(f"max_age ({self.max_age}) must be positive")
        if self.max_age > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"max_age ({self.max_age} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma)"
            )
        if not 0 < self.tail_prob < 1:
            raise ValueError(f"tail_prob ({self.tail_prob}) must be between 0 and 1")

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        elif target_software in ["r8s", "pyr8s"]:
            return f"constrain taxon={{label}} max_age={self.max_age};"
        elif target_software == "treepl":
            return f"max = {{Name}} {self.max_age}"
        elif target_software == "lsd2":
            # LSD2 date 文件取值轴：value = -age（负值表示距今），``v1:v2`` 表示
            # value ∈ [v1, v2]（NA 表示该侧无界）。经 IQ-TREE 2.3.6 / 3.1.3 实测
            # 回读校准节点年龄确认：``NA:-T`` 被 LSD 解码为“年龄 ≥ T”（下界），
            # ``-T:NA`` 解码为“年龄 ≤ T”（上界）。
            # 上界约束：年龄至多 max_age → 下端（更年轻侧）开放 -> -max_age:NA
            return f"-{self.max_age}:NA"
        elif target_software == "pathd8":
            # PATHd8 支持单独 maxage
            return f"mrca: {{t1}}, {{t2}}, maxage={self.max_age};"
        else:
            # 审阅项 C-1：旧实现在这里先警告"降级为 [min_val, max_age] 均匀区间"，
            # 紧接着对 wlogdate/mdcat 返回 max_age/2 的**点估计**，两条警告互相矛盾、
            # 且第一条与返回值不符（min_val 只用于拼字符串）。现在每个分支只发一条
            # 与实际返回值一致的警告；未知目标软件不再发警告而是直接报错。
            if target_software in ["wlogdate", "mdcat"]:
                # wLogDate/MD-Cat 仅支持点估计。为与 UniformAgeConstraint 保持统一
                # 的降级原则（区间中点），MaximumAgeConstraint 视为 [0, max_age]
                # 的退化区间，取中点 max_age/2 作为点估计。
                mid = self.max_age / 2.0
                _emit_degradation_warning(
                    f"{target_software} only supports point calibrations and does not "
                    f"support maximum-only constraints. Degrading the "
                    f"≤{self.max_age} Ma upper bound to the point estimate {mid} Ma "
                    f"(midpoint of the implied range [0.0, {self.max_age}] Ma). The "
                    f"upper bound itself is NOT enforced: {target_software} may return "
                    f"an age older than {self.max_age} Ma."
                )
                return f"{{Name}}={{taxa_str}}\t{mid}"
            else:
                raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """MCMCTree ``U(tU, tailR)``（单位 Ga）。

        上游 UPPER_F 只写 2 个槽（treesub.c:8740-8743 → pfossil[1]=b、
        pfossil[2]=tailR），与本串逐位吻合（审阅报告 §六.2）。
        """
        return _render_mcmctree_upper_bound(
            source=f"MaximumAgeConstraint(max_age={self.max_age} Ma)",
            upper_ma=self.max_age,
            tail_upper=self.tail_prob,
        )

    def describe(self) -> str:
        """人类可读摘要：最大年龄（上界）约束。"""
        return f"≤{self.max_age} Ma"


@dataclass
class SoftBoundsConstraint(AgeConstraint):
    """双软边界约束（MCMCTree ``B(tL, tU, tailL, tailR)`` 格式）"""

    min_age: float  # 最小年龄（Ma）
    max_age: float  # 最大年龄（Ma）
    # 审阅项 B-1：BOUND_F 的第 3/4 槽分别是 tailL/tailR，必须分开表达；旧实现只有
    # 一个 tail_prob，写进第 3 槽后 tailR 悄悄停在默认 0.025，"对称软界"变成非对称。
    tail_lower: float = MCMCTREE_DEFAULT_TAIL_PROB  # 区间左尾概率（槽 3）
    tail_upper: float = MCMCTREE_DEFAULT_TAIL_PROB  # 区间右尾概率（槽 4）
    # 旧字段/便捷入口：对称软界。给出时等价于同时设 tail_lower=tail_upper=tail_prob。
    tail_prob: Optional[float] = None

    def __post_init__(self) -> None:
        _require_finite_number(self.min_age, "min_age")
        _require_finite_number(self.max_age, "max_age")
        if self.min_age <= 0:
            raise ValueError(f"min_age ({self.min_age}) must be positive")
        if self.max_age > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"max_age ({self.max_age} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma)"
            )
        if self.min_age >= self.max_age:
            raise ValueError(
                f"min_age ({self.min_age}) must be < max_age ({self.max_age})"
            )
        if self.tail_prob is not None:
            _require_finite_number(self.tail_prob, "tail_prob")
            if self.tail_lower == MCMCTREE_DEFAULT_TAIL_PROB and self.tail_upper == (
                MCMCTREE_DEFAULT_TAIL_PROB
            ):
                # 只给了旧的对称字段 tail_prob → 对称地填进两个槽
                self.tail_lower = self.tail_prob
                self.tail_upper = self.tail_prob
            elif self.tail_lower != self.tail_prob or self.tail_upper != self.tail_prob:
                raise ValueError(
                    "SoftBoundsConstraint: tail_prob (symmetric convenience field) "
                    f"({self.tail_prob}) conflicts with the explicit "
                    f"tail_lower={self.tail_lower}/tail_upper={self.tail_upper}; use "
                    "tail_lower/tail_upper only"
                )
        _require_finite_number(self.tail_lower, "tail_lower")
        _require_finite_number(self.tail_upper, "tail_upper")
        for name, value in (
            ("tail_lower", self.tail_lower),
            ("tail_upper", self.tail_upper),
        ):
            if not 0 < value < 1:
                raise ValueError(
                    f"{name} ({value}) must be a probability strictly between 0 and 1"
                )
        if self.tail_lower + self.tail_upper >= 1:
            raise ValueError(
                f"tail_lower + tail_upper ({self.tail_lower + self.tail_upper}) "
                f"must be < 1"
            )

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        else:
            # 降级为均匀分布（wlogdate/mdcat 只能落点估计，见下方分支的专用警告）
            if target_software not in ("wlogdate", "mdcat"):
                _emit_degradation_warning(
                    f"{target_software} does not support soft bounds. "
                    f"Degrading to uniform constraint [{self.min_age}, "
                    f"{self.max_age}] (tail probabilities tailL={self.tail_lower}, "
                    f"tailR={self.tail_upper} are lost)"
                )
            if target_software in ["r8s", "pyr8s"]:
                return (
                    f"constrain taxon={{label}} min_age={self.min_age};\n"
                    f"constrain taxon={{label}} max_age={self.max_age};"
                )
            elif target_software == "treepl":
                return (
                    f"min = {{Name}} {self.min_age}\n" f"max = {{Name}} {self.max_age}"
                )
            elif target_software == "pathd8":
                return (
                    f"mrca: {{t1}}, {{t2}}, minage={self.min_age};\n"
                    f"mrca: {{t1}}, {{t2}}, maxage={self.max_age};"
                )
            elif target_software == "lsd2":
                return f"-{self.max_age}:-{self.min_age}"
            elif target_software in ["wlogdate", "mdcat"]:
                mid = (self.min_age + self.max_age) / 2
                _emit_degradation_warning(
                    f"{target_software} only supports point calibrations. Degrading "
                    f"the soft bounds [{self.min_age}, {self.max_age}] Ma to the "
                    f"midpoint {mid}; both bounds and the tail probabilities "
                    f"(tailL={self.tail_lower}, tailR={self.tail_upper}) are lost."
                )
                return f"{{Name}}={{taxa_str}}\t{mid}"
            else:
                raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """MCMCTree ``B(tL, tU, tailL, tailR)``（单位 Ga，4 个槽写满，审阅项 B-1）。"""
        return _render_mcmctree_bounds(
            source=(
                f"SoftBoundsConstraint(min_age={self.min_age} Ma, "
                f"max_age={self.max_age} Ma)"
            ),
            lower_ma=self.min_age,
            upper_ma=self.max_age,
            tail_lower=self.tail_lower,
            tail_upper=self.tail_upper,
        )

    def describe(self) -> str:
        """人类可读摘要：双软边界约束。"""
        return (
            f"{self.min_age}-{self.max_age} Ma (soft, "
            f"tailL={self.tail_lower}, tailR={self.tail_upper})"
        )


@dataclass
class GammaPriorConstraint(AgeConstraint):
    """Gamma 先验约束（MCMCTree 的 ``G(alpha, beta)``，单位 Ga）

    **上游只接受 2 个参数**（审阅项 A-4）：``npfossils[GAMMA_F] = 2``
    （mcmctree.c:202-204），解析式为
    `sscanf(pch, "%lf,%lf", &pfossil[0], &pfossil[1])`（treesub.c:8754-8756，4.8 与
    4.10.8 完全相同），密度 `lnpt = a*log(b) - b*t + (a-1)*log(t) - lgamma(a)`
    （mcmctree.c:2622-2624）**没有平移项**、支撑从 0 开始。因此本类的 ``offset``
    在 MCMCTree 通路上无法表达：写出第三个值连被读入的机会都没有。渲染时会发
    :class:`SemanticDegradationWarning` 明确披露这一点。
    """

    alpha: float  # 形状参数
    beta: float  # 速率参数（1/Ma）
    offset: float = 0.0  # 偏移量（Ma）——MCMCTree 通路不生效，见类 docstring
    scale: float = 1.0  # 尺度参数（MCMCTree 无此参数，通过并入速率 beta 表达）

    def __post_init__(self) -> None:
        _require_finite_number(self.alpha, "alpha")
        _require_finite_number(self.beta, "beta")
        _require_finite_number(self.scale, "scale")
        _require_finite_number(self.offset, "offset")
        if self.alpha <= 0:
            raise ValueError(f"alpha ({self.alpha}) must be positive")
        if self.beta <= 0:
            raise ValueError(f"beta ({self.beta}) must be positive")
        if self.scale <= 0:
            raise ValueError(f"scale ({self.scale}) must be positive")
        if self.offset < 0:
            raise ValueError(f"offset ({self.offset}) must be non-negative")
        if self.offset > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"offset ({self.offset} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma - Earth's age)"
            )

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        else:
            # 降级为均匀分布：从 Gamma 分布计算均值和 95% 区间
            # MCMCTree 参数化: X = offset + scale * Gamma(alpha, beta)
            # scipy.stats.gamma(shape=k, scale=θ) 的 ppf 返回 scale=θ 尺度的分位数
            # Gamma(alpha, beta) 的 scipy 参数化: shape=alpha, scale=1/beta
            # 因此 X = offset + scale * Gamma(alpha,beta) 的分位数:
            #   ppf_X(q) = offset + scale * stats.gamma.ppf(q, alpha, scale=1/beta)
            # 注意：这条通路上 offset 是真的参与计算的（与 MCMCTree 通路不同）。
            scipy_scale = 1.0 / self.beta
            ci_lower = (
                stats.gamma.ppf(0.025, self.alpha, scale=scipy_scale) * self.scale
                + self.offset
            )
            ci_upper = (
                stats.gamma.ppf(0.975, self.alpha, scale=scipy_scale) * self.scale
                + self.offset
            )
            ci_lower, ci_upper, notes = _validate_derived_interval(
                f"GammaPriorConstraint(alpha={self.alpha}, beta={self.beta})",
                ci_lower,
                ci_upper,
            )

            _emit_degradation_warning(
                f"{target_software} does not support Gamma priors. "
                f"Degrading to uniform constraint [{ci_lower:.1f}, {ci_upper:.1f}] "
                f"(derived from Gamma({self.alpha}, {self.beta}, {self.offset}, scale={self.scale}))"
                + _notes_suffix(notes)
            )
            if target_software in ["r8s", "pyr8s"]:
                return (
                    f"constrain taxon={{label}} min_age={ci_lower:.1f};\n"
                    f"constrain taxon={{label}} max_age={ci_upper:.1f};"
                )
            elif target_software == "treepl":
                return (
                    f"min = {{Name}} {ci_lower:.1f}\n" f"max = {{Name}} {ci_upper:.1f}"
                )
            elif target_software == "pathd8":
                return (
                    f"mrca: {{t1}}, {{t2}}, minage={ci_lower:.1f};\n"
                    f"mrca: {{t1}}, {{t2}}, maxage={ci_upper:.1f};"
                )
            elif target_software == "lsd2":
                return f"-{ci_upper:.1f}:-{ci_lower:.1f}"
            elif target_software in ["wlogdate", "mdcat"]:
                mean = self.alpha * self.scale / self.beta + self.offset
                return f"{{Name}}={{taxa_str}}\t{mean:.1f}"
            else:
                raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """MCMCTree ``G(alpha, beta_Ga)``（**仅 2 个参数**，审阅项 A-4）。

        上游事实（4.8 与 4.10.8 一致）：

        * ``npfossils[GAMMA_F] = 2``（mcmctree.c:202-204）；
        * ``sscanf(pch, "%lf,%lf", &pfossil[0], &pfossil[1])``（treesub.c:8754-8756）
          ——第三个值**连被读入的机会都没有**；
        * 密度 ``lnpt = a*log(b) - b*t + (a-1)*log(t) - lgamma(a)``
          （mcmctree.c:2622-2624）是标准 ``Gamma(α, rate=β)``，**没有平移项**，
          支撑从 0 开始。

        因此旧串 ``G(alpha, beta, offset/1000)`` 里的 offset 是被静默丢弃的第 3 参数，
        旧注释"MCMCTree 仅支持 3 参数 G(alpha, beta, offset)"与上游矛盾。现在：

        1. 只写 2 个参数，不再声称存在可用的 offset 槽；
        2. ``offset != 0`` 时发 :class:`SemanticDegradationWarning`，说明"化石地板"
           在 MCMCTree 通路上没有生效，并给出可用的替代表达（``B()``/``U()``/
           ``SoftBoundsConstraint`` 以 offset 为下界）；
        3. 量纲换算本身保持不变（审阅报告 §六.4 确认正确）：密度以 Ga 为 t 单位，
           ``E[t_Ma] = alpha/beta`` 换算到 Ga 需要把速率乘 1000；
           ``c·Gamma(α,β) = Gamma(α, β/c)``，故 ``scale`` 通过 ``beta*1000/scale``
           并入速率，先验均值/方差与 ``offset + scale·Gamma(alpha, beta)``（offset=0
           时）完全一致。
        """
        effective_beta = self.beta * 1000.0 / self.scale
        if not math.isfinite(effective_beta) or effective_beta <= 0:
            raise ValueError(
                f"GammaPriorConstraint(alpha={self.alpha}, beta={self.beta}, "
                f"scale={self.scale}): rendered MCMCTree rate parameter "
                f"({effective_beta}) must be a positive finite number"
            )
        if self.offset != 0:
            _emit_degradation_warning(
                f"MCMCTree parses only 2 parameters for G() "
                f"(npfossils[GAMMA_F] = 2, treesub.c:8754-8756) and its Gamma density "
                f"has no shift term (mcmctree.c:2622-2624), so the offset/floor of "
                f"{self.offset} Ma is NOT applied: the rendered prior "
                f"G({self.alpha}, {effective_beta:.6f}) has support starting at 0 and "
                f"its mode can be arbitrarily close to 0. To express a hard/soft "
                f"stratigraphic floor, use SoftBoundsConstraint (or "
                f"UniformAgeConstraint) with min_age={self.offset} Ma instead of a "
                f"shifted Gamma."
            )
        return f"G({self.alpha}, {effective_beta:.6f})"

    def describe(self) -> str:
        """人类可读摘要：Gamma 先验约束（标注 offset 在 MCMCTree 通路未生效）。"""
        offset_note = (
            ", NOT applied by MCMCTree" if self.offset else ", no stratigraphic floor"
        )
        return (
            f"Gamma(α={self.alpha}, β={self.beta}, "
            f"offset={self.offset} Ma{offset_note}, scale={self.scale})"
        )


@dataclass
class SkewNormalConstraint(AgeConstraint):
    """偏态正态先验（MCMCTree 4.8+）"""

    location: float  # 位置参数（Ma）
    scale: float  # 尺度参数（Ma）
    shape: float  # 形状参数（偏度）

    def __post_init__(self) -> None:
        _require_finite_number(self.location, "location")
        _require_finite_number(self.scale, "scale")
        _require_finite_number(self.shape, "shape")
        if self.scale <= 0:
            raise ValueError(f"scale ({self.scale}) must be positive")
        if self.location < 0:
            raise ValueError(f"location ({self.location}) must be non-negative")
        if self.location > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"location ({self.location} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma - Earth's age)"
            )

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        else:
            # 降级为均匀分布：使用 scipy.stats 计算实际分位数
            method_notes = []
            try:
                from scipy.stats import skewnorm

                ci_lower = skewnorm.ppf(
                    0.025, self.shape, loc=self.location, scale=self.scale
                )
                ci_upper = skewnorm.ppf(
                    0.975, self.shape, loc=self.location, scale=self.scale
                )
                mean = skewnorm.mean(self.shape, loc=self.location, scale=self.scale)
            except Exception:
                # 回退：使用矩估计 + 对称区间近似
                method_notes.append(
                    "scipy.stats.skewnorm was unavailable, so the interval is a "
                    "symmetric normal approximation (mean ± 1.96·sd) and not the true "
                    "equal-tailed 95% interval of a skewed prior"
                )
                delta = self.shape / (1 + self.shape**2) ** 0.5
                mean = self.location + self.scale * delta * (2 / math.pi) ** 0.5
                std = self.scale * (1 - 2 * delta**2 / math.pi) ** 0.5
                ci_lower = mean - 1.96 * std
                ci_upper = mean + 1.96 * std

            # 审阅项 C-5：负下界的 0 截断改变区间含义，必须写进警告文本
            ci_lower, ci_upper, interval_notes = _validate_derived_interval(
                f"SkewNormalConstraint(location={self.location}, "
                f"scale={self.scale}, shape={self.shape})",
                ci_lower,
                ci_upper,
            )
            notes = interval_notes + method_notes
            if target_software in ("wlogdate", "mdcat"):
                notes = [
                    f"{target_software} only supports point calibrations, so the "
                    f"prior mean is written and the interval is dropped entirely"
                ] + notes

            _emit_degradation_warning(
                f"{target_software} does not support Skew Normal priors. "
                + (
                    f"Writing the prior mean {float(mean):.1f} Ma."
                    if target_software in ("wlogdate", "mdcat")
                    else f"Degrading to uniform constraint "
                    f"[{ci_lower:.1f}, {ci_upper:.1f}]"
                )
                + f" (derived from SN({self.location}, {self.scale}, {self.shape}))"
                + _notes_suffix(notes)
            )
            if target_software in ["r8s", "pyr8s"]:
                return (
                    f"constrain taxon={{label}} min_age={ci_lower:.1f};\n"
                    f"constrain taxon={{label}} max_age={ci_upper:.1f};"
                )
            elif target_software == "treepl":
                return (
                    f"min = {{Name}} {ci_lower:.1f}\n" f"max = {{Name}} {ci_upper:.1f}"
                )
            elif target_software == "pathd8":
                return (
                    f"mrca: {{t1}}, {{t2}}, minage={ci_lower:.1f};\n"
                    f"mrca: {{t1}}, {{t2}}, maxage={ci_upper:.1f};"
                )
            elif target_software == "lsd2":
                return f"-{ci_upper:.1f}:-{ci_lower:.1f}"
            elif target_software in ["wlogdate", "mdcat"]:
                return f"{{Name}}={{taxa_str}}\t{mean:.1f}"
            else:
                raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """MCMCTree ``SN(location, scale, shape)``（单位 Ga，3 个槽）。

        上游 SKEWN_F 恰好读 3 个值（``npfossils[SKEWN_F] = 3``，
        treesub.c:8757-8759 → ``logPDFSkewN(t, p[0], p[1], p[2])``），与本串对应
        （审阅报告 §六.3 确认）。
        """
        source = f"SkewNormalConstraint(location={self.location}, scale={self.scale})"
        return (
            f"SN({_render_ga_bound(self.location, 'location', source)}, "
            f"{_render_ga_bound(self.scale, 'scale', source)}, "
            f"{_format_number(self.shape)})"
        )

    def describe(self) -> str:
        """人类可读摘要：偏态正态先验约束。"""
        return (
            f"SkewNormal(loc={self.location}, scale={self.scale}, "
            f"shape={self.shape})"
        )


@dataclass
class SkewTConstraint(AgeConstraint):
    """偏态 t 分布先验（MCMCTree 4.8+）"""

    location: float  # 位置参数（Ma）
    scale: float  # 尺度参数（Ma）
    shape: float  # 形状参数（偏度）
    df: float  # 自由度

    def __post_init__(self) -> None:
        _require_finite_number(self.location, "location")
        _require_finite_number(self.scale, "scale")
        _require_finite_number(self.shape, "shape")
        _require_finite_number(self.df, "df")
        if self.scale <= 0:
            raise ValueError(f"scale ({self.scale}) must be positive")
        if self.df <= 2:
            raise ValueError(
                f"df ({self.df}) must be > 2 for t-distribution to have finite variance"
            )
        if self.location < 0:
            raise ValueError(f"location ({self.location}) must be non-negative")
        if self.location > MAX_GEOLOGICAL_AGE_MA:
            raise ValueError(
                f"location ({self.location} Ma) exceeds geological reasonable bound "
                f"({MAX_GEOLOGICAL_AGE_MA} Ma - Earth's age)"
            )

    def moment_approximation(self) -> Tuple[float, float]:
        """Azzalini skew-t 的**精确**前两阶矩 + 正态近似 95% 区间（审阅项 C-3/C-4）。

        .. math::

            \\delta = \\lambda / \\sqrt{1+\\lambda^2},\\quad
            E[X] = \\xi + \\omega\\,\\delta \\sqrt{\\nu/\\pi}\\,
                   \\Gamma((\\nu-1)/2)/\\Gamma(\\nu/2)

            Var[X] = \\omega^2 \\left[ \\frac{\\nu}{\\nu-2}(1+\\delta^2)
                     - \\mu_{shift}^2 \\right], \\qquad \\nu > 2

        说明（C-3）：本方法**不**尝试提供精确分位数。scipy 从未提供 skew-t 分布
        （``scipy.stats.skewt`` 不存在；``skewnorm``/``gamma``/``t``/``nct`` 存在），
        旧代码 ``from scipy.stats import skewt`` 因此恒抛 ``ImportError``，"精确分位数"
        分支永不可达、又从不告知用户。可选出路是改用 ``scipy.stats.nct``（非中心 t）
        并实现 Azzalini skew-union 表示，但那条映射并非逐参数等价、且与 MCMCTree 自己
        ``PDFSkewT``（mcmctree.c:2629）所用定义是否同构需要上游数值核对——风险高于
        收益，因此这里显式放弃"精确"分支，只保留有据可依的矩公式，并把近似方式与
        误差方向写进降级警告。方差公式（C-4）不再有 ``df <= 2`` 的死分支：
        ``__post_init__`` 已经要求 ``df > 2``，所以有限方差恒成立。
        """
        delta = self.shape / (1 + self.shape**2) ** 0.5
        try:
            gamma_ratio = math.gamma((self.df - 1) / 2) / math.gamma(self.df / 2)
        except (ValueError, OverflowError):  # 仅当 df 大到使 math.gamma 溢出时可达
            gamma_ratio = 1.0
        mean_shift = delta * (self.df / math.pi) ** 0.5 * gamma_ratio
        mean = self.location + self.scale * mean_shift
        var = self.scale**2 * (self.df / (self.df - 2) * (1 + delta**2) - mean_shift**2)
        if not math.isfinite(var) or var < 0:
            raise ValueError(
                f"SkewTConstraint(location={self.location}, scale={self.scale}, "
                f"shape={self.shape}, df={self.df}): variance is not a finite "
                f"non-negative number ({var})"
            )
        std = var**0.5
        return mean, std

    def _render_software_format(self, target_software: str) -> str:
        if target_software == "mcmctree":
            return self.to_mcmctree_calib_string()
        else:
            # 降级为均匀分布：使用 Azzalini skew-t 的精确均值/方差 + 正态近似区间
            mean, std = self.moment_approximation()
            ci_lower = mean - 1.96 * std
            ci_upper = mean + 1.96 * std

            # 审阅项 C-5：0 截断改变区间含义（不再是等尾 95%），必须记入警告文本
            ci_lower, ci_upper, interval_notes = _validate_derived_interval(
                f"SkewTConstraint(location={self.location}, scale={self.scale}, "
                f"shape={self.shape}, df={self.df})",
                ci_lower,
                ci_upper,
            )
            notes = interval_notes + [
                "the interval is a symmetric NORMAL approximation "
                "mean ± 1.96·sd built from the exact skew-t mean/variance, because "
                "scipy provides no skew-t distribution (scipy.stats.skewt does not "
                "exist) so the true equal-tailed quantiles cannot be computed; the "
                "real skew-t is asymmetric and heavier-tailed (excess kurtosis "
                "6/(df-4) for df>4), so this interval is too narrow on the heavy-tail "
                f"side — the error grows as df approaches 2 (now df={self.df})"
            ]
            if target_software in ("wlogdate", "mdcat"):
                notes = [
                    f"{target_software} only supports point calibrations, so the "
                    f"prior mean is written and the interval is dropped entirely"
                ] + notes

            _emit_degradation_warning(
                f"{target_software} does not support Skew T priors. "
                + (
                    f"Writing the prior mean {mean:.1f} Ma."
                    if target_software in ("wlogdate", "mdcat")
                    else f"Degrading to uniform constraint "
                    f"[{ci_lower:.1f}, {ci_upper:.1f}]"
                )
                + " (derived from "
                f"ST({self.location}, {self.scale}, {self.shape}, {self.df}))"
                + _notes_suffix(notes)
            )
            if target_software in ["r8s", "pyr8s"]:
                return (
                    f"constrain taxon={{label}} min_age={ci_lower:.1f};\n"
                    f"constrain taxon={{label}} max_age={ci_upper:.1f};"
                )
            elif target_software == "treepl":
                return (
                    f"min = {{Name}} {ci_lower:.1f}\n" f"max = {{Name}} {ci_upper:.1f}"
                )
            elif target_software == "pathd8":
                return (
                    f"mrca: {{t1}}, {{t2}}, minage={ci_lower:.1f};\n"
                    f"mrca: {{t1}}, {{t2}}, maxage={ci_upper:.1f};"
                )
            elif target_software == "lsd2":
                return f"-{ci_upper:.1f}:-{ci_lower:.1f}"
            elif target_software in ["wlogdate", "mdcat"]:
                return f"{{Name}}={{taxa_str}}\t{mean:.1f}"
            else:
                raise ValueError(f"Unknown software: {target_software}")

    def to_mcmctree_calib_string(self) -> str:
        """MCMCTree ``ST(location, scale, shape, df)``（单位 Ga，4 个槽）。

        上游 SKEWT_F 恰好读 4 个值（``npfossils[SKEWT_F] = 4``，
        treesub.c:8760-8762 → ``log(PDFSkewT(t, p[0], p[1], p[2], p[3]))``），
        与本串逐位吻合（审阅报告 §六.3 确认）。
        """
        source = f"SkewTConstraint(location={self.location}, scale={self.scale})"
        return (
            f"ST({_render_ga_bound(self.location, 'location', source)}, "
            f"{_render_ga_bound(self.scale, 'scale', source)}, "
            f"{_format_number(self.shape)}, {_format_number(self.df)})"
        )

    def describe(self) -> str:
        """人类可读摘要：偏态 t 分布先验约束。"""
        return (
            f"SkewT(loc={self.location}, scale={self.scale}, "
            f"shape={self.shape}, df={self.df})"
        )
