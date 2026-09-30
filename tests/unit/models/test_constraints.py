"""
年龄约束类型的单元测试

测试所有约束类的创建、验证和格式转换功能
"""

import inspect
import math
import re
import warnings

import pytest

from phylodater.core.exceptions import SemanticDegradationWarning
from phylodater.models.calibration import CalibrationPoint, FossilMetadata
from phylodater.models.constraints import (
    MAX_GEOLOGICAL_AGE_MA,
    FixedAgeConstraint,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    SkewNormalConstraint,
    SkewTConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)


class TestFixedAgeConstraint:
    """测试固定年龄约束"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = FixedAgeConstraint(fixed_age=50.0)
        assert constraint.fixed_age == 50.0

    def test_integer_age(self):
        """测试整数年龄"""
        constraint = FixedAgeConstraint(fixed_age=100)
        assert constraint.fixed_age == 100

    def test_zero_age_raises(self):
        """测试零年龄应抛出异常"""
        with pytest.raises(ValueError, match="must be positive"):
            FixedAgeConstraint(fixed_age=0)

    def test_negative_age_raises(self):
        """测试负年龄应抛出异常"""
        with pytest.raises(ValueError, match="must be positive"):
            FixedAgeConstraint(fixed_age=-10)

    def test_excessive_age_raises(self):
        """测试超过地球年龄应抛出异常"""
        with pytest.raises(ValueError, match="exceeds geological reasonable bound"):
            FixedAgeConstraint(fixed_age=MAX_GEOLOGICAL_AGE_MA + 1)

    def test_non_numeric_age_raises(self):
        """测试非数字年龄应抛出异常"""
        with pytest.raises(TypeError, match="must be a number"):
            FixedAgeConstraint(fixed_age="50")

    def test_nan_age_raises(self):
        """NaN 必须被拒（审阅项 B-19：旧实现里 NaN 与任何数比较均为 False 而全绿通过）"""
        with pytest.raises(ValueError, match="must be a finite real number"):
            FixedAgeConstraint(fixed_age=float("nan"))

    def test_infinite_age_raises(self):
        """±inf 必须被拒（与 B-19 同一道守卫；旧实现靠上限比较偶然拦住 +inf）"""
        with pytest.raises(ValueError, match="must be a finite real number"):
            FixedAgeConstraint(fixed_age=float("inf"))
        with pytest.raises(ValueError, match="must be a finite real number"):
            FixedAgeConstraint(fixed_age=float("-inf"))

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式（审阅项 A-3/B-1：4 个槽 + 渲染后严格下界<上界）"""
        constraint = FixedAgeConstraint(fixed_age=1000.0)
        result = constraint.to_mcmctree_calib_string()
        # MCMCtree 不使用点校准（@X 格式），而是用窄的软边界近似
        # 1000.0 Ma = 1.0 Ga，容差 1% 即 0.01 Ga；尾部概率每侧 0.005（合 99% 在界内）
        assert result == "B(0.9900, 1.0100, 0.005, 0.005)"

    def test_to_mcmctree_format_renders_four_slot_bounds(self):
        """B() 必须写满 4 个值（审阅项 B-1：第 3 槽是 tailL，tailR 不能被默认值顶替）"""
        parts = FixedAgeConstraint(fixed_age=100.0).to_mcmctree_calib_string()[2:-1]
        assert len(parts.split(",")) == 4

    def test_to_r8s_format(self):
        """测试转换为 r8s 格式"""
        constraint = FixedAgeConstraint(fixed_age=50.0)
        result = constraint.to_software_format("r8s")
        assert "fixage" in result
        assert "50.0" in result

    def test_to_treepl_format(self):
        """测试转换为 treePL 格式"""
        constraint = FixedAgeConstraint(fixed_age=50.0)
        result = constraint.to_software_format("treepl")
        assert "min" in result
        assert "max" in result
        assert "50.0" in result

    def test_to_lsd2_format(self):
        """测试转换为 LSD2 格式"""
        constraint = FixedAgeConstraint(fixed_age=50.0)
        result = constraint.to_software_format("lsd2")
        assert result == "-50.0"

    def test_to_unknown_software_raises(self):
        """测试未知软件应抛出异常"""
        constraint = FixedAgeConstraint(fixed_age=50.0)
        with pytest.raises(ValueError, match="Unknown software"):
            constraint.to_software_format("unknown")


class TestUniformAgeConstraint:
    """测试均匀分布年龄约束"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = UniformAgeConstraint(min_age=10.0, max_age=50.0)
        assert constraint.min_age == 10.0
        assert constraint.max_age == 50.0

    def test_equal_min_max_raises(self):
        """测试 min >= max 应抛出异常"""
        with pytest.raises(ValueError, match="must be < max_age"):
            UniformAgeConstraint(min_age=50.0, max_age=50.0)

    def test_min_greater_than_max_raises(self):
        """测试 min > max 应抛出异常"""
        with pytest.raises(ValueError, match="must be < max_age"):
            UniformAgeConstraint(min_age=60.0, max_age=50.0)

    def test_negative_min_raises(self):
        """测试负 min_age 应抛出异常"""
        with pytest.raises(ValueError, match="must be positive"):
            UniformAgeConstraint(min_age=-10.0, max_age=50.0)

    def test_excessive_max_raises(self):
        """测试超过地球年龄的 max 应抛出异常"""
        with pytest.raises(ValueError, match="exceeds geological reasonable bound"):
            UniformAgeConstraint(min_age=10.0, max_age=MAX_GEOLOGICAL_AGE_MA + 1)

    def test_non_numeric_raises(self):
        """测试非数字应抛出异常"""
        with pytest.raises(TypeError, match="must be a number"):
            UniformAgeConstraint(min_age="10", max_age=50.0)

    def test_nan_bounds_raise(self):
        """测试 NaN 界被拒（审阅项 B-19）"""
        with pytest.raises(ValueError, match="must be a finite real number"):
            UniformAgeConstraint(min_age=float("nan"), max_age=100.0)
        with pytest.raises(ValueError, match="must be a finite real number"):
            UniformAgeConstraint(min_age=10.0, max_age=float("nan"))

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式（审阅项 B-1：4 槽，且不再是会除零的 p=0）"""
        constraint = UniformAgeConstraint(min_age=500.0, max_age=1000.0)
        result = constraint.to_mcmctree_calib_string()
        assert result == "B(0.5000, 1.0000, 0.0001, 0.0001)"

    def test_to_lsd2_format(self):
        """测试转换为 LSD2 格式"""
        constraint = UniformAgeConstraint(min_age=10.0, max_age=50.0)
        result = constraint.to_software_format("lsd2")
        assert result == "-50.0:-10.0"

    def test_to_wlogdate_degradation(self):
        """测试 wLogDate 降级警告"""
        constraint = UniformAgeConstraint(min_age=10.0, max_age=50.0)
        with pytest.warns(SemanticDegradationWarning):
            result = constraint.to_software_format("wlogdate")
        assert result is not None
        assert "30.0" in result  # midpoint


class TestSoftLowerBoundConstraint:
    """测试软下界约束（审阅项 B-2：字段与 L() 的 4 个槽一一对应）"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = SoftLowerBoundConstraint(min_age=50.0)
        assert constraint.min_age == 50.0
        # 槽 2 = 相对偏移比例 P、槽 3 = Cauchy 尺度 c、槽 4 = 真正的尾概率
        assert constraint.offset_fraction == 0.1
        assert constraint.cauchy_scale == 1.0
        assert constraint.tail_prob == 0.025

    def test_custom_shape_params(self):
        """测试自定义截断 Cauchy 参数"""
        constraint = SoftLowerBoundConstraint(
            min_age=50.0, offset_fraction=0.05, cauchy_scale=0.5, tail_prob=0.1
        )
        assert constraint.offset_fraction == 0.05
        assert constraint.cauchy_scale == 0.5
        assert constraint.tail_prob == 0.1

    def test_legacy_tail_prob_tail_shape_still_maps(self):
        """旧 API (tail_prob=P, tail_shape=c) 仍可用，数值语义不变（DeprecationWarning）"""
        with pytest.warns(DeprecationWarning, match="deprecated"):
            constraint = SoftLowerBoundConstraint(
                min_age=50.0, tail_prob=0.05, tail_shape=0.5
            )
        # 旧的 tail_prob 其实是上游的 P，旧的 tail_shape 其实是上游的 c
        assert constraint.offset_fraction == 0.05
        assert constraint.cauchy_scale == 0.5
        assert constraint.tail_prob == 0.025
        assert constraint.to_mcmctree_calib_string() == "L(0.0500, 0.05, 0.5, 0.025)"

    def test_invalid_tail_prob_raises(self):
        """测试无效尾概率应抛出异常（槽 4 必须是 (0,1) 开区间）"""
        with pytest.raises(ValueError, match="must be a probability strictly"):
            SoftLowerBoundConstraint(min_age=50.0, tail_prob=0)
        with pytest.raises(ValueError, match="must be a probability strictly"):
            SoftLowerBoundConstraint(min_age=50.0, tail_prob=1.0)

    def test_invalid_cauchy_scale_raises(self):
        """测试非正尺度应抛出异常"""
        with pytest.raises(ValueError, match="cauchy_scale"):
            SoftLowerBoundConstraint(min_age=50.0, cauchy_scale=0)
        with pytest.raises(ValueError, match="offset_fraction"):
            SoftLowerBoundConstraint(min_age=50.0, offset_fraction=-0.1)

    def test_nan_age_raises(self):
        """测试 NaN 被拒（审阅项 B-19）"""
        with pytest.raises(ValueError, match="must be a finite real number"):
            SoftLowerBoundConstraint(min_age=float("nan"))

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式：L(tL, P, c, tailL) 4 个槽写满（审阅项 B-2）"""
        constraint = SoftLowerBoundConstraint(
            min_age=500.0, offset_fraction=0.1, cauchy_scale=1.0, tail_prob=0.025
        )
        result = constraint.to_mcmctree_calib_string()
        assert result == "L(0.5000, 0.1, 1.0, 0.025)"
        # 第 4 槽（真正的尾概率）过去从不写出，永远停在上游默认 0.025
        assert len(result[2:-1].split(",")) == 4

    def test_to_lsd2_format(self):
        """测试转换为 LSD2 格式（软下界 -> NA:-T，经 IQ-TREE 2.x/3.x 实测）"""
        constraint = SoftLowerBoundConstraint(min_age=50.0)
        result = constraint.to_software_format("lsd2")
        assert result == "NA:-50.0"

    def test_to_other_software_degradation(self):
        """测试其他软件降级警告"""
        constraint = SoftLowerBoundConstraint(min_age=50.0)
        with pytest.warns(SemanticDegradationWarning):
            result = constraint.to_software_format("r8s")
        assert result is not None
        assert "50" in result


class TestMaximumAgeConstraint:
    """测试最大年龄约束"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = MaximumAgeConstraint(max_age=100.0)
        assert constraint.max_age == 100.0
        assert constraint.tail_prob == 0.025

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式（U 只有 2 个槽，与上游一致）"""
        constraint = MaximumAgeConstraint(max_age=1000.0, tail_prob=0.025)
        result = constraint.to_mcmctree_calib_string()
        assert result == "U(1.0000, 0.025)"

    def test_nan_age_raises(self):
        """测试 NaN 被拒（审阅项 B-19）"""
        with pytest.raises(ValueError, match="must be a finite real number"):
            MaximumAgeConstraint(max_age=float("nan"))

    def test_wlogdate_degradation_warning_matches_returned_value(self):
        """审阅项 C-1：降级分支只能有一条与实际返回值一致的警告

        旧实现先警告"降级为 [min_val, max_age] 均匀区间"，随后返回 max_age/2 的点
        估计，两条警告互相矛盾且第一条与返回值不符。
        """
        constraint = MaximumAgeConstraint(max_age=100.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = constraint.to_software_format("wlogdate", name="C", taxa_str="A+B")
        messages = [str(w.message) for w in caught]
        assert len(messages) == 1, f"应只有一条降级警告，实际 {messages}"
        assert "point estimate 50.0" in messages[0]
        assert "uniform range" not in messages[0]
        assert result == "C=A+B\t50.0"
        # 不再留下只用于拼字符串的 min_val：警告文本里不得出现 [10.0, 100.0] 之类的假区间
        assert "[10.0, 100.0]" not in messages[0]

    def test_unknown_software_raises_without_warning(self):
        """未知目标软件直接报错，不再先发一条与结果无关的警告"""
        constraint = MaximumAgeConstraint(max_age=100.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with pytest.raises(ValueError, match="Unknown software"):
                constraint.to_software_format("no-such-tool", name="C", taxa_str="A+B")
        assert [w for w in caught if w.category is SemanticDegradationWarning] == []

    def test_to_lsd2_format(self):
        """测试转换为 LSD2 格式（上界 -> -T:NA，经 IQ-TREE 2.x/3.x 实测）"""
        constraint = MaximumAgeConstraint(max_age=100.0)
        result = constraint.to_software_format("lsd2")
        assert result == "-100.0:NA"


class TestSoftBoundsConstraint:
    """测试双软边界约束（审阅项 B-1：tailL/tailR 必须分开写满 4 个槽）"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = SoftBoundsConstraint(min_age=10.0, max_age=100.0)
        assert constraint.min_age == 10.0
        assert constraint.max_age == 100.0
        # 两侧尾部各自独立，默认值与上游一致（treesub.c:8675）
        assert constraint.tail_lower == 0.025
        assert constraint.tail_upper == 0.025

    def test_symmetric_tail_prob_alias(self):
        """旧字段 tail_prob 作为"对称软界"入口仍可用，且两侧同时生效"""
        constraint = SoftBoundsConstraint(min_age=10.0, max_age=100.0, tail_prob=0.1)
        assert constraint.tail_lower == 0.1
        assert constraint.tail_upper == 0.1
        assert constraint.to_mcmctree_calib_string() == "B(0.0100, 0.1000, 0.1, 0.1)"

    def test_asymmetric_tails_are_honoured(self):
        """审阅项 B-1：旧的单个 tail_prob 会让 tailR 悄悄停在 0.025，现在两侧独立"""
        constraint = SoftBoundsConstraint(
            min_age=500.0, max_age=1000.0, tail_lower=0.1, tail_upper=0.01
        )
        assert constraint.to_mcmctree_calib_string() == "B(0.5000, 1.0000, 0.1, 0.01)"

    def test_conflicting_tail_fields_raise(self):
        """tail_prob 与不一致的 tail_lower/tail_upper 同时给出时报错而不是静默取舍"""
        with pytest.raises(ValueError, match="conflicts"):
            SoftBoundsConstraint(
                min_age=10.0, max_age=100.0, tail_lower=0.05, tail_prob=0.1
            )

    def test_tails_must_sum_below_one(self):
        with pytest.raises(ValueError, match="must be < 1"):
            SoftBoundsConstraint(
                min_age=10.0, max_age=100.0, tail_lower=0.6, tail_upper=0.6
            )

    def test_nan_age_raises(self):
        """测试 NaN 被拒（审阅项 B-19）"""
        with pytest.raises(ValueError, match="must be a finite real number"):
            SoftBoundsConstraint(min_age=10.0, max_age=float("nan"))

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式"""
        constraint = SoftBoundsConstraint(
            min_age=500.0, max_age=1000.0, tail_prob=0.025
        )
        result = constraint.to_mcmctree_calib_string()
        assert result == "B(0.5000, 1.0000, 0.025, 0.025)"

    def test_to_other_software_degradation(self):
        """测试其他软件降级警告"""
        constraint = SoftBoundsConstraint(min_age=10.0, max_age=100.0)
        with pytest.warns(SemanticDegradationWarning):
            result = constraint.to_software_format("r8s")
        assert result is not None
        assert "10" in result and "100" in result


class TestGammaPriorConstraint:
    """测试 Gamma 先验约束"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = GammaPriorConstraint(alpha=2.0, beta=0.5)
        assert constraint.alpha == 2.0
        assert constraint.beta == 0.5
        assert constraint.offset == 0.0
        assert constraint.scale == 1.0

    def test_invalid_alpha_raises(self):
        """测试无效 alpha 应抛出异常"""
        with pytest.raises(ValueError, match="must be positive"):
            GammaPriorConstraint(alpha=0, beta=0.5)

    def test_invalid_beta_raises(self):
        """测试无效 beta 应抛出异常"""
        with pytest.raises(ValueError, match="must be positive"):
            GammaPriorConstraint(alpha=2.0, beta=0)

    def test_nan_params_raise(self):
        """测试 NaN 被拒（审阅项 B-19）"""
        with pytest.raises(ValueError, match="must be a finite real number"):
            GammaPriorConstraint(alpha=float("nan"), beta=0.5)
        with pytest.raises(ValueError, match="must be a finite real number"):
            GammaPriorConstraint(alpha=2.0, beta=0.5, offset=float("nan"))

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式（审阅项 A-4：G() 只有 2 个参数，offset 不被上游读取）

        旧断言 ``"G(2.0, 500.000000, 0.1000)"`` 编码的正是被报告判为错误的行为：
        ``npfossils[GAMMA_F] = 2``，``sscanf`` 只读两个值，密度里也没有平移项，
        所以第 3 个值连被读入的机会都没有。
        """
        constraint = GammaPriorConstraint(alpha=2.0, beta=0.5, offset=100.0)
        with pytest.warns(SemanticDegradationWarning, match="is NOT applied"):
            result = constraint.to_mcmctree_calib_string()
        assert result == "G(2.0, 500.000000)"
        assert len(result[2:-1].split(",")) == 2

    def test_to_mcmctree_format_without_offset_does_not_warn(self):
        """offset == 0 时不该出现"地板丢失"的警告；速率换算保持不变"""
        constraint = GammaPriorConstraint(alpha=2.0, beta=0.5, scale=2.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = constraint.to_mcmctree_calib_string()
        assert [w for w in caught if w.category is SemanticDegradationWarning] == []
        # effective_beta = beta*1000/scale = 0.5*1000/2 = 250
        assert result == "G(2.0, 250.000000)"

    def test_describe_discloses_inactive_offset(self):
        """describe() 必须披露 offset 在 MCMCTree 通路未生效（审阅项 A-4）"""
        desc = GammaPriorConstraint(alpha=2.0, beta=0.5, offset=34.0).describe()
        assert "NOT applied by MCMCTree" in desc
        assert "34.0" in desc

    def test_offset_still_used_by_non_mcmctree_degradation(self):
        """降级到区间时 offset 仍然参与计算（只有 MCMCTree 通路无法表达平移）"""
        constraint = GammaPriorConstraint(alpha=2.0, beta=0.5, offset=100.0)
        with pytest.warns(SemanticDegradationWarning):
            result = constraint.to_software_format("r8s", label="T")
        assert "min_age=" in result
        # 95% 区间下界 = 100 + Gamma(2,0.5) 的 2.5% 分位数 > 100
        lower = float(re.search(r"min_age=([0-9.]+)", result).group(1))
        assert lower > 100.0

    def test_to_other_software_degradation(self):
        """测试其他软件降级警告"""
        constraint = GammaPriorConstraint(alpha=2.0, beta=0.5)
        with pytest.warns(SemanticDegradationWarning):
            result = constraint.to_software_format("r8s")
        assert result is not None


class TestSkewNormalConstraint:
    """测试偏态正态约束"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = SkewNormalConstraint(location=50.0, scale=10.0, shape=2.0)
        assert constraint.location == 50.0
        assert constraint.scale == 10.0
        assert constraint.shape == 2.0

    def test_invalid_scale_raises(self):
        """测试无效 scale 应抛出异常"""
        with pytest.raises(ValueError, match="must be positive"):
            SkewNormalConstraint(location=50.0, scale=0, shape=2.0)

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式"""
        constraint = SkewNormalConstraint(location=500.0, scale=100.0, shape=2.0)
        result = constraint.to_mcmctree_calib_string()
        assert result == "SN(0.5000, 0.1000, 2.0)"


class TestSkewTConstraint:
    """测试偏态 t 分布约束"""

    def test_valid_creation(self):
        """测试有效值创建"""
        constraint = SkewTConstraint(location=50.0, scale=10.0, shape=2.0, df=5.0)
        assert constraint.location == 50.0
        assert constraint.scale == 10.0
        assert constraint.shape == 2.0
        assert constraint.df == 5.0

    def test_invalid_df_raises(self):
        """测试无效自由度应抛出异常"""
        with pytest.raises(ValueError, match="must be > 2"):
            SkewTConstraint(location=50.0, scale=10.0, shape=2.0, df=2.0)

    def test_to_mcmctree_format(self):
        """测试转换为 MCMCTree 格式"""
        constraint = SkewTConstraint(location=500.0, scale=100.0, shape=2.0, df=5.0)
        result = constraint.to_mcmctree_calib_string()
        assert result == "ST(0.5000, 0.1000, 2.0, 5.0)"

    # ---- 审阅项 C-3 / C-4 / C-5 ------------------------------------------- #

    def test_no_fake_precise_quantile_branch(self):
        """scipy 从来没有 skewt：不得再有"看似精确、实则永不可达"的分位数分支"""
        import scipy.stats

        assert not hasattr(scipy.stats, "skewt")
        source = inspect.getsource(SkewTConstraint)
        assert "import skewt as" not in source
        assert "skewt_dist" not in source

    def test_degradation_warning_states_approximation_and_error(self):
        """降级警告必须写清正态近似方式与误差方向（C-3：过去静默回退）"""
        constraint = SkewTConstraint(location=100.0, scale=10.0, shape=5.0, df=3.5)
        with pytest.warns(SemanticDegradationWarning) as caught:
            constraint.to_software_format("r8s", label="T")
        msg = str(caught[0].message)
        assert "NORMAL approximation" in msg
        assert "scipy.stats.skewt does not exist" in msg
        assert "too narrow on the heavy-tail side" in msg

    def test_moment_approximation_matches_documented_formulas(self):
        """矩公式与 Azzalini skew-t 一致，且方差含 (1+δ²) 因子（C-4：无死分支）"""
        c = SkewTConstraint(location=100.0, scale=10.0, shape=5.0, df=10.0)
        delta = c.shape / (1 + c.shape**2) ** 0.5
        gamma_ratio = math.gamma((c.df - 1) / 2) / math.gamma(c.df / 2)
        mean_shift = delta * (c.df / math.pi) ** 0.5 * gamma_ratio
        mean, std = c.moment_approximation()
        assert mean == pytest.approx(c.location + c.scale * mean_shift)
        var = c.scale**2 * (c.df / (c.df - 2) * (1 + delta**2) - mean_shift**2)
        assert std == pytest.approx(var**0.5)
        # 方差恒按 df>2 的公式计算（__post_init__ 已排除 df<=2），无 std=scale 回退
        assert std != c.scale
        source = inspect.getsource(SkewTConstraint)
        assert "std = self.scale" not in source

    def test_truncated_lower_bound_is_disclosed(self):
        """C-5：下界被 0 截断时必须写进警告文本，且不能再自称等尾 95% 区间"""
        # location 小而尺度大 → 2.5% 分位数为负
        constraint = SkewTConstraint(location=5.0, scale=10.0, shape=0.0, df=3.0)
        with pytest.warns(SemanticDegradationWarning) as caught:
            result = constraint.to_software_format("treepl", name="C")
        msg = str(caught[0].message)
        assert "truncated" in msg and "equal-tailed" in msg
        lower = float(re.search(r"min = \S+ ([0-9.]+)", result).group(1))
        assert lower == 0.0

    def test_skew_normal_truncation_is_disclosed(self):
        """C-5 同族：SN 的负下界截断同样要留痕"""
        constraint = SkewNormalConstraint(location=1.0, scale=10.0, shape=-5.0)
        with pytest.warns(SemanticDegradationWarning) as caught:
            constraint.to_software_format("r8s", label="T")
        msg = str(caught[0].message)
        assert "truncated" in msg

    def test_skew_normal_wlogdate_warning_matches_point_output(self):
        """wlogdate/mdcat 落点估计时，警告文本不得声称写出了区间（C-1 同族）"""
        constraint = SkewNormalConstraint(location=50.0, scale=10.0, shape=2.0)
        with pytest.warns(SemanticDegradationWarning) as caught:
            result = constraint.to_software_format("wlogdate", name="C", taxa_str="A+B")
        msg = str(caught[0].message)
        assert "only supports point calibrations" in msg
        assert result.startswith("C=A+B\t")
        mean = float(result.split("\t")[1])
        assert f"mean {mean:.1f} Ma" in msg


class TestRenderedIntervalGuards:
    """审阅项 A-3 / C-2：格式化精度与"最终渲染值"守卫"""

    @pytest.mark.parametrize("fixed_age", [0.005, 0.02, 0.05, 0.1, 0.5, 1.0, 10.0])
    def test_fixed_age_rendered_bounds_strictly_ordered(self, fixed_age):
        """报告 A-3 表中的 7 个年龄：打印后必须严格 lower < upper"""
        rendered = FixedAgeConstraint(fixed_age=fixed_age).to_mcmctree_calib_string()
        parts = rendered[2:-1].split(",")
        assert len(parts) == 4
        lower, upper = float(parts[0]), float(parts[1])
        assert upper > lower, f"{fixed_age} Ma 渲染后区间退化: {rendered}"
        assert lower > 0.0, f"{fixed_age} Ma 渲染后下界非正: {rendered}"
        # 两个界都关于固定年龄对称（±1%）
        assert (lower + upper) / 2 == pytest.approx(fixed_age / 1000.0, rel=1e-6)

    def test_fixed_age_tenth_of_a_ma_no_longer_inverted(self):
        """旧实现输出 B(0.0001, 0.0000)：下界 > 上界，MCMCTree 直接中止"""
        assert (
            FixedAgeConstraint(fixed_age=0.1).to_mcmctree_calib_string()
            == "B(0.000099, 0.000101, 0.005, 0.005)"
        )

    @pytest.mark.parametrize("pair", [(0.05, 0.1), (6.0, 8.0), (1e-4, 2e-4)])
    def test_uniform_rendered_bounds_strictly_ordered(self, pair):
        lo, hi = pair
        rendered = UniformAgeConstraint(
            min_age=lo, max_age=hi
        ).to_mcmctree_calib_string()
        parts = rendered[2:-1].split(",")
        assert float(parts[1]) > float(parts[0])

    @pytest.mark.parametrize("fixed_age", [4560.0, 4600.0])
    def test_rendered_upper_bound_cannot_exceed_earth_age(self, fixed_age):
        """C-2：输入 4600 Ma 通过校验，但 ±1% 窗口把上界推到 4.646 Ga > 地球年龄"""
        with pytest.raises(ValueError, match="exceeds the maximum geological age"):
            FixedAgeConstraint(fixed_age=fixed_age).to_mcmctree_calib_string()

    def test_rendered_interval_error_for_degenerate_case(self):
        """构造一个渲染后坍缩的极端区间：必须报错而不是产出 B(a, a, ...)"""
        constraint = UniformAgeConstraint(min_age=1e-12, max_age=1.0000001e-12)
        with pytest.raises(ValueError, match="degenerate"):
            constraint.to_mcmctree_calib_string()

    def test_nan_cannot_reach_rendered_string(self):
        """B-19 的出口侧兜底：即使绕过建构期校验也不得产出含 nan 的串"""
        constraint = UniformAgeConstraint(min_age=10.0, max_age=100.0)
        object.__setattr__(constraint, "min_age", float("nan"))
        with pytest.raises(ValueError, match="not a finite number"):
            constraint.to_mcmctree_calib_string()


class TestCalibrationPointBounds:
    """审阅项 C-6 / C-7：CalibrationPoint 与 FossilMetadata 的边界与交叉校验"""

    @pytest.mark.parametrize(
        "constraint, expected_lower, expected_upper",
        [
            (FixedAgeConstraint(fixed_age=100.0), True, True),
            (UniformAgeConstraint(min_age=50.0, max_age=100.0), True, True),
            (SoftBoundsConstraint(min_age=50.0, max_age=100.0), True, True),
            (SoftLowerBoundConstraint(min_age=50.0), True, False),
            (MaximumAgeConstraint(max_age=100.0), False, True),
        ],
    )
    def test_bound_flags(self, constraint, expected_lower, expected_upper):
        point = CalibrationPoint(name="n", age_constraint=constraint)
        assert point.has_lower_bound is expected_lower
        # C-7：固定点同时是上界与下界
        assert point.has_upper_bound is expected_upper

    def test_fixed_point_counts_as_upper_bound(self):
        point = CalibrationPoint(
            name="n", age_constraint=FixedAgeConstraint(fixed_age=100.0)
        )
        assert point.has_upper_bound and point.has_lower_bound

    def test_fossil_metadata_rejects_reversed_stratigraphic_range(self):
        """C-6：min > max 过去被静默接受"""
        with pytest.raises(ValueError, match="must be <="):
            FossilMetadata(min_stratigraphic_age=120.0, max_stratigraphic_age=100.0)
        with pytest.raises(ValueError, match="must be <="):
            FossilMetadata.from_dict(
                {"min_stratigraphic_age": 120, "max_stratigraphic_age": 100}
            )

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), -1.0, 5000.0])
    def test_fossil_metadata_rejects_invalid_ages(self, bad):
        with pytest.raises(ValueError):
            FossilMetadata.from_dict({"min_stratigraphic_age": bad})

    def test_fossil_metadata_allows_empty_and_ordered_range(self):
        assert FossilMetadata().min_stratigraphic_age is None
        meta = FossilMetadata(min_stratigraphic_age=100.0, max_stratigraphic_age=100.0)
        assert meta.max_stratigraphic_age == 100.0

    def test_disjoint_fossil_and_constraint_warns(self):
        """C-6：地层区间与校准区间毫不重叠时必须有警告（方向约定留给作者核实）"""
        with pytest.warns(UserWarning, match="do not overlap at all"):
            CalibrationPoint(
                name="n",
                age_constraint=FixedAgeConstraint(fixed_age=50.0),
                fossil_metadata=FossilMetadata(
                    min_stratigraphic_age=100.0, max_stratigraphic_age=120.0
                ),
            )

    def test_overlapping_fossil_and_constraint_is_silent(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            CalibrationPoint(
                name="n",
                age_constraint=UniformAgeConstraint(min_age=150.0, max_age=200.0),
                fossil_metadata=FossilMetadata(
                    min_stratigraphic_age=100.0, max_stratigraphic_age=120.0
                ),
            )
        assert [w for w in caught if "do not overlap" in str(w.message)] == []

    def test_from_dict_propagates_metadata_validation(self):
        payload = CalibrationPoint(
            name="n", age_constraint=UniformAgeConstraint(min_age=10.0, max_age=20.0)
        ).to_dict()
        payload["fossil_metadata"] = {
            "min_stratigraphic_age": 15.0,
            "max_stratigraphic_age": 14.0,
        }
        with pytest.raises(ValueError, match="must be <="):
            CalibrationPoint.from_dict(payload)
