"""
回归测试：针对 code-review/phylodater-code-review-2026-07-23.md 中已修复的关键 bug。

本文件聚焦「纯单元测试」，不依赖任何外部定年二进制
（mcmctree / r8s / treePL / lsd2 / pathd8 等均未安装，相关集成/功能测试应 skip/xfail）。

覆盖项（与任务清单一一对应）：
  - constraints.SkewTConstraint 方差公式含 (1+δ²)
  - constraints.GammaPriorConstraint.to_mcmctree_calib_string 输出 2 参数 G(...)
    （offset 上游无槽可读，改为警告/describe 披露；审阅项 A-4）
  - constraints.UniformAgeConstraint.to_mcmctree_calib_string 输出 4 参数
    B(tL, tU, tailL, tailR)，尾部极小但非零（审阅项 B-1）
  - constraints 全部 __post_init__ 拒绝 NaN 年龄（审阅项 B-19）
  - 固定年龄的 4600 Ma 守卫作用在**最终渲染值**上（审阅项 C-2）
  - 8 个约束子类均有可用的 describe()
  - tree.get_rooting_info() 有根/无根判定 + 含 likely_artifact/shortest_ratio
  - tree.is_monophyletic() 子集匹配（不误用子串包含）
  - tree.get_mrca() 找不到时返回 None
  - tree.get_distance_matrix() 缓存生效且结果一致
  - tree.write(format="nexus") 产出 NEXUS 标记
  - calibration_resolver.resolve_special_nodes 将 Methanobacteriota 归入 bacteria
  - calibration_resolver._is_ancestor 在 MRCA 为 None/异常时返回 False
  - auto_calibrator.get_calibration_summary 对 Gamma/Skew 不再显示 unknown
  - mcmctree_method RootAge 单位已修正为 Ga（+ DEFAULT_ROOT_AGE_GA 常量）
"""

import math
import re
import warnings
from types import SimpleNamespace

import pytest

from phylodater.models.calibration import CalibrationPoint
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
from phylodater.models.tree import PhylogeneticTree
from phylodater.services.auto_calibrator import AutoCalibrator
from phylodater.services.calibration_resolver import (
    CalibrationResolver,
    resolve_special_nodes,
)

# 8 个约束子类（有序），用于 describe() 通测
_ALL_CONSTRAINTS = [
    FixedAgeConstraint(fixed_age=100.0),
    UniformAgeConstraint(min_age=50.0, max_age=100.0),
    SoftLowerBoundConstraint(min_age=50.0),
    MaximumAgeConstraint(max_age=100.0),
    SoftBoundsConstraint(min_age=50.0, max_age=100.0),
    GammaPriorConstraint(alpha=2.0, beta=0.5, offset=10.0, scale=1.0),
    SkewNormalConstraint(location=100.0, scale=10.0, shape=0.5),
    SkewTConstraint(location=100.0, scale=10.0, shape=0.5, df=10.0),
]


# --------------------------------------------------------------------------- #
# 1) SkewTConstraint 方差公式必须含 (1+δ²) 因子
# --------------------------------------------------------------------------- #
class TestSkewTConstraintVariance:
    def test_variance_formula_includes_one_plus_delta_squared(self, monkeypatch):
        """内联方差公式必须含 (1+δ²)，与漏该因子的旧值明显不同。

        scipy 1.17 的 scipy.stats 不含 skewt，to_software_format 走 fallback 分支
        并触发内联方差计算；这里强制将 skewt 置为 None 以稳定复现该路径
        （与 scipy 版本无关）。
        """
        import scipy.stats as sp_stats

        monkeypatch.setattr(sp_stats, "skewt", None, raising=False)

        c = SkewTConstraint(location=100.0, scale=10.0, shape=5.0, df=10.0)

        delta = c.shape / (1 + c.shape**2) ** 0.5
        gamma_ratio = math.gamma((c.df - 1) / 2) / math.gamma(c.df / 2)
        mean_shift = delta * (c.df / math.pi) ** 0.5 * gamma_ratio
        mean = c.location + c.scale * mean_shift

        # 修复后的方差公式：scale²·[ν/(ν-2)·(1+δ²) − mean_shift²]
        var_fixed = c.scale**2 * (c.df / (c.df - 2) * (1 + delta**2) - mean_shift**2)
        # 旧（有 bug）公式：漏掉 (1+δ²)
        var_buggy = c.scale**2 * (c.df / (c.df - 2) - mean_shift**2)

        std_fixed = var_fixed**0.5
        std_buggy = var_buggy**0.5

        ci_lo_fixed = max(0.0, mean - 1.96 * std_fixed)
        ci_hi_fixed = mean + 1.96 * std_fixed
        ci_lo_buggy = max(0.0, mean - 1.96 * std_buggy)
        ci_hi_buggy = mean + 1.96 * std_buggy

        _expected_fixed = (
            f"constrain taxon={{label}} min_age={ci_lo_fixed:.1f};\n"
            f"constrain taxon={{label}} max_age={ci_hi_fixed:.1f}"
        )
        expected_buggy = (
            f"constrain taxon={{label}} min_age={ci_lo_buggy:.1f};\n"
            f"constrain taxon={{label}} max_age={ci_hi_buggy:.1f}"
        )

        actual = c.to_software_format("r8s")

        # 1) 解析实际输出中的 ci 边界，断言其数值等于「修复后」公式推导结果
        lo = float(re.search(r"min_age=([0-9.]+)", actual).group(1))
        hi = float(re.search(r"max_age=([0-9.]+)", actual).group(1))
        assert lo == pytest.approx(ci_lo_fixed, abs=0.05)
        assert hi == pytest.approx(ci_hi_fixed, abs=0.05)
        # 2) 与漏 (1+δ²) 的旧结果不同 → 证明修复生效
        assert actual != expected_buggy
        # 3) 数值上方差必须显著包含 (1+δ²) 因子（df=10, shape=5 时约为 3×）
        assert var_fixed > 1.5 * var_buggy


# --------------------------------------------------------------------------- #
# 2) GammaPriorConstraint.to_mcmctree_calib_string → 2 参数 G(...)（审阅项 A-4）
# --------------------------------------------------------------------------- #
class TestGammaPriorConstraintMcmcTree:
    def test_two_param_g_with_effective_beta_and_offset_disclosure(self):
        """G(alpha, effective_beta)：**只有 2 个参数**，offset 不上串但必须留痕。

        原断言"输出 3 参数 G(alpha, beta, offset)"编码的正是审阅报告判为错误的行为：
        上游 ``npfossils[GAMMA_F] = 2``、``sscanf(pch, "%lf,%lf", ...)``
        （treesub.c:8754-8756，4.8 与 4.10.8 相同），Gamma 密度里也没有平移项
        （mcmctree.c:2622-2624），所以第 3 个值连被读入的机会都没有。

        effective_beta = beta·1000 / scale 的量纲换算本身是**正确**的（报告 §六.4
        已确认），保持不变；offset 改为在警告与 describe() 中披露"未生效"。
        """
        from phylodater.core.exceptions import SemanticDegradationWarning

        c = GammaPriorConstraint(alpha=2.0, beta=4.0, offset=500.0, scale=2.0)
        with pytest.warns(SemanticDegradationWarning, match="is NOT applied"):
            s = c.to_mcmctree_calib_string()

        assert s.startswith("G(")
        assert s.endswith(")")
        # 2 个参数 => 恰好 1 个逗号（第 3 个 offset 槽在上游不存在）
        assert s.count(",") == 1, f"expected 2 params, got: {s}"

        alpha_s, beta_s = s[2:-1].split(",")
        assert float(alpha_s) == pytest.approx(2.0)
        # effective_beta = beta·1000/scale = 4.0·1000/2.0 = 2000.0
        assert float(beta_s) == pytest.approx(2000.0)
        # offset 不再被写成第 3 参数，但必须在人类可读摘要里披露"未生效"
        assert "500.0" in c.describe()
        assert "NOT applied by MCMCTree" in c.describe()

    def test_zero_offset_renders_without_warning(self):
        from phylodater.core.exceptions import SemanticDegradationWarning

        c = GammaPriorConstraint(alpha=2.0, beta=4.0, scale=2.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            s = c.to_mcmctree_calib_string()
        assert s == "G(2.0, 2000.000000)"
        assert [w for w in caught if w.category is SemanticDegradationWarning] == []


# ----------------------------------------------------------------------- #
# 3) UniformAgeConstraint.to_mcmctree_calib_string → 4 参数 B(...)（审阅项 B-1）
# --------------------------------------------------------------------------- #
class TestUniformAgeConstraintMcmcTree:
    def test_four_param_b_with_non_zero_tails(self):
        """B(tL, tU, tailL, tailR)：4 个槽写满，且第 3 槽不再是会除零的 0.0

        原断言"3 参数 B(..., 0.0)"编码的是被报告判为错误的行为：BOUND_F 期望 4 个
        值（``npfossils[BOUND_F] = 4``），第 3 槽是 **tailL**、tailR 停在默认 0.025；
        传 0 会让 ``thetaL = (1-tailL-tailR)*a/(tailL*(b-a))`` 除零、
        ``log(tailL*thetaL/a)`` 成为 NaN，"硬界"右侧还悄悄留着 2.5% 的尾巴。
        """
        c = UniformAgeConstraint(min_age=500.0, max_age=1000.0)
        s = c.to_mcmctree_calib_string()

        assert s.startswith("B(")
        parts = s[2:-1].split(",")
        # 4 个参数 => 3 个逗号
        assert len(parts) == 4, f"expected 4 params, got: {s}"
        assert s.count(",") == 3
        # min/max 已 ÷1000
        assert float(parts[0]) == pytest.approx(0.5)
        assert float(parts[1]) == pytest.approx(1.0)
        # 两个尾部概率都是极小但**非零**的值（近似硬界，且不触发除零）
        for tail in parts[2:]:
            value = float(tail)
            assert 0.0 < value < 0.5, f"尾部概率非法: {s}"
        assert float(parts[1]) > float(parts[0])

    def test_soft_bounds_tails_are_separate_slots(self):
        """SoftBoundsConstraint 的 tailL/tailR 各自独立（审阅项 B-1 第 2 条）"""
        c = SoftBoundsConstraint(
            min_age=500.0, max_age=1000.0, tail_lower=0.1, tail_upper=0.025
        )
        assert c.to_mcmctree_calib_string() == "B(0.5000, 1.0000, 0.1, 0.025)"

    def test_soft_lower_bound_emits_four_slot_L(self):
        """L(tL, P, c, tailL)：第 4 槽（真正的尾概率）现在写得出门（审阅项 B-2）"""
        c = SoftLowerBoundConstraint(
            min_age=500.0, offset_fraction=0.1, cauchy_scale=1.0, tail_prob=0.05
        )
        s = c.to_mcmctree_calib_string()
        assert s == "L(0.5000, 0.1, 1.0, 0.05)"

    def test_nan_ages_are_rejected_at_construction(self):
        """审阅项 B-19：NaN 过去穿过全部 >/< 校验，产出 B(0.0001, nan, 0.01)"""
        for cls, kwargs in [
            (FixedAgeConstraint, {"fixed_age": math.nan}),
            (UniformAgeConstraint, {"min_age": math.nan, "max_age": 100.0}),
            (UniformAgeConstraint, {"min_age": 10.0, "max_age": math.nan}),
            (SoftBoundsConstraint, {"min_age": 10.0, "max_age": math.nan}),
            (SoftLowerBoundConstraint, {"min_age": math.nan}),
            (MaximumAgeConstraint, {"max_age": math.nan}),
            (GammaPriorConstraint, {"alpha": math.nan, "beta": 1.0}),
            (SkewNormalConstraint, {"location": 10.0, "scale": math.nan, "shape": 1.0}),
            (
                SkewTConstraint,
                {"location": 10.0, "scale": 1.0, "shape": math.nan, "df": 5.0},
            ),
        ]:
            with pytest.raises(ValueError, match="finite"):
                cls(**kwargs)

    def test_fixed_age_guard_validates_rendered_value(self):
        """审阅项 C-2：输入 4600 Ma 合法，但 ±1% 窗口的上界 4.646 Ga 必须被拦下"""
        assert FixedAgeConstraint(fixed_age=4554.0).to_mcmctree_calib_string()
        with pytest.raises(ValueError, match="maximum geological age"):
            FixedAgeConstraint(
                fixed_age=MAX_GEOLOGICAL_AGE_MA
            ).to_mcmctree_calib_string()


# --------------------------------------------------------------------------- #
# 4) 8 个约束子类都有 describe() 且返回有效字符串
# --------------------------------------------------------------------------- #
class TestAllConstraintsHaveDescribe:
    @pytest.mark.parametrize(
        "constraint", _ALL_CONSTRAINTS, ids=lambda c: type(c).__name__
    )
    def test_describe_returns_nonempty_string(self, constraint):
        desc = constraint.describe()
        assert isinstance(desc, str)
        assert desc.strip() != ""

    def test_describe_is_distinct_per_constraint_type(self):
        descs = [c.describe() for c in _ALL_CONSTRAINTS]
        # 8 种约束应产生 8 个不同描述（多态生效，不会全部回退到同一字符串）
        assert len(set(descs)) == len(_ALL_CONSTRAINTS)

    def test_base_describe_contains_type_name(self):
        # 基类默认实现应返回类名（验证 describe 机制本身含「类型名」）
        from phylodater.models.constraints import AgeConstraint

        class _Dummy(AgeConstraint):
            def _render_software_format(self, target_software):
                return ""

            def to_mcmctree_calib_string(self):
                return ""

        assert _Dummy().describe() == "_Dummy"


# --------------------------------------------------------------------------- #
# 5) tree.get_rooting_info()：有根/无根判定 + 含必要键
# --------------------------------------------------------------------------- #
class TestTreeRooting:
    def test_rooted_tree_reports_rooted_with_required_keys(self):
        # 根节点二叉 -> 已定根
        t = PhylogeneticTree("((A,B),(C,D));")
        info = t.get_rooting_info()

        assert info["is_rooted"] is True
        # 必须包含任务要求的键
        for key in ("is_rooted", "likely_artifact", "shortest_ratio"):
            assert key in info
        assert isinstance(info["shortest_ratio"], float)
        assert isinstance(info["likely_artifact"], bool)

    def test_unrooted_tree_reports_unrooted(self):
        # 根节点三叉（multifurcation）-> 未定根
        t = PhylogeneticTree("((A,B),(C,D),E);")
        info = t.get_rooting_info()

        assert info["is_rooted"] is False
        assert "likely_artifact" in info
        assert "shortest_ratio" in info


# --------------------------------------------------------------------------- #
# 6) tree.is_monophyletic()：子集匹配，不误用子串包含
# --------------------------------------------------------------------------- #
class TestTreeMonophyly:
    def test_is_monophyletic_no_substring_false_match(self):
        # 含 "A" 与 "AB" 两个独立叶；二者是姐妹群，MRCA 恰好为 {A, AB}
        t = PhylogeneticTree("((A,AB),C);")

        # "A" 与 "AB" 按精确叶名集合构成单系（不为子串 "A"⊂"AB" 所混淆）
        assert t.is_monophyletic(["A", "AB"]) is True
        # 真实拓扑：{"A","C"} 的 MRCA 含 {A,AB,C} ≠ {"A","C"} → 非单系
        assert t.is_monophyletic(["A", "C"]) is False
        # 单叶在 ETE3 严格语义下非单系；返回 bool 且不抛错（兜底路径稳定）
        assert isinstance(t.is_monophyletic(["A"]), bool)


# --------------------------------------------------------------------------- #
# 7) tree.get_mrca()：找不到时返回 None（而非整个叶列表）
# --------------------------------------------------------------------------- #
class TestTreeMrca:
    def test_get_mrca_returns_none_when_not_found(self):
        t = PhylogeneticTree("((A,B),(C,D));")
        # 不存在的 tip 应使 get_mrca 返回 None
        assert t.get_mrca(["Z", "Y"]) is None

    def test_get_mrca_returns_leaf_names_when_found(self):
        t = PhylogeneticTree("((A,B),(C,D));")
        mrca = t.get_mrca(["A", "B"])
        assert set(mrca) == {"A", "B"}


# --------------------------------------------------------------------------- #
# 8) tree.get_distance_matrix()：缓存生效且结果一致
# --------------------------------------------------------------------------- #
class TestTreeDistanceMatrix:
    def test_distance_matrix_cached_and_consistent(self):
        t = PhylogeneticTree("((A:1,B:1):1,(C:1,D:1):1);")
        m1 = t.get_distance_matrix()
        m2 = t.get_distance_matrix()

        # 第二次调用应返回同一缓存对象（缓存生效）
        assert m2 is m1
        assert len(m1) > 0
        # 对称性
        assert m1[("A", "B")] == m1[("B", "A")]
        # A-B 距离为 2（共享父枝 1+1）
        assert m1[("A", "B")] == pytest.approx(2.0)


# --------------------------------------------------------------------------- #
# 9) tree.write(format="nexus")：产出 NEXUS 标记
# --------------------------------------------------------------------------- #
class TestTreeWriteNexus:
    def test_write_nexus_contains_nexus_markers(self, tmp_path):
        t = PhylogeneticTree("((A:1,B:1),(C:1,D:1));")
        out = tmp_path / "tree.nex"
        t.write(out, format="nexus")

        content = out.read_text(encoding="utf-8")
        assert "begin trees;" in content.lower()
        assert "tree" in content.lower()


# --------------------------------------------------------------------------- #
# 10) calibration_resolver.resolve_special_nodes：Methanobacteriota 归入 bacteria
# --------------------------------------------------------------------------- #
class TestResolveSpecialNodes:
    def test_methanobacteriota_classified_as_bacteria(self):
        newick = (
            "((Bacteria_phylum_Methanobacteriota_sp1,Bacteria_phylum_Other_sp2),"
            "(Archaea_sp3,Archaea_sp4));"
        )
        t = PhylogeneticTree(newick)
        resolver = CalibrationResolver(t)  # 无 taxonomy_parser

        result = resolve_special_nodes(t, resolver)

        assert "LBCA" in result
        assert "LUCA" in result
        assert "LACA" in result

        # Methanobacteriota 必须归入 bacteria（LBCA/LUCA），而非被错误排除
        methano = "Bacteria_phylum_Methanobacteriota_sp1"
        assert methano in result["LBCA"]
        assert methano in result["LUCA"]
        # 不应出现在古菌列表
        assert methano not in result["LACA"]

    def test_methanobacteria_substring_not_excluded(self):
        # 代码注释所述场景：含 "methanobacteria" 子串的 tip 不应被错误排除
        newick = "((Methanobacteria_abc,Eubacteria_def)," "(Archaea_ghi,Archaea_jkl));"
        t = PhylogeneticTree(newick)
        resolver = CalibrationResolver(t)
        result = resolve_special_nodes(t, resolver)

        assert "Methanobacteria_abc" in result["LBCA"]
        assert "Methanobacteria_abc" in result["LUCA"]


# --------------------------------------------------------------------------- #
# 11) calibration_resolver._is_ancestor：MRCA 为 None/异常时返回 False
# --------------------------------------------------------------------------- #
class TestIsAncestorNone:
    def test_is_ancestor_false_when_mrca_none(self):
        t = PhylogeneticTree("((A,B),(C,D));")
        resolver = CalibrationResolver(t)

        # 幽灵分类单元（不在树中）→ get_mrca 返回 None → 应返回 False（非 True）
        cal1 = CalibrationPoint(
            name="g1", resolved_taxa=["GhostX", "GhostY"], is_root_node=False
        )
        cal2 = CalibrationPoint(name="g2", resolved_taxa=["GhostX"], is_root_node=False)
        assert resolver._is_ancestor(cal1, cal2) is False

    def test_is_ancestor_true_for_real_ancestry(self):
        t = PhylogeneticTree("(A,(B,C));")
        resolver = CalibrationResolver(t)

        cal1 = CalibrationPoint(
            name="rootish", resolved_taxa=["A", "B", "C"], is_root_node=False
        )
        cal2 = CalibrationPoint(name="bc", resolved_taxa=["B", "C"], is_root_node=False)
        assert resolver._is_ancestor(cal1, cal2) is True

    def test_is_ancestor_false_when_not_ancestor(self):
        t = PhylogeneticTree("(A,(B,C));")
        resolver = CalibrationResolver(t)

        cal1 = CalibrationPoint(name="bc", resolved_taxa=["B", "C"], is_root_node=False)
        cal2 = CalibrationPoint(name="ab", resolved_taxa=["A", "B"], is_root_node=False)
        assert resolver._is_ancestor(cal1, cal2) is False


# --------------------------------------------------------------------------- #
# 12) auto_calibrator.get_calibration_summary：Gamma/Skew 不再显示 unknown
# --------------------------------------------------------------------------- #
class TestAutoCalibratorSummary:
    def test_summary_uses_describe_for_gamma_and_skew(self):
        t = PhylogeneticTree("((A,B),(C,D));")
        ac = AutoCalibrator(t)

        cal_gamma = CalibrationPoint(
            name="g",
            age_constraint=GammaPriorConstraint(
                alpha=2.0, beta=0.5, offset=10.0, scale=1.0
            ),
            resolved_taxa=["A", "B"],
            mrca_leaf_pair=("A", "B"),
        )
        cal_skew = CalibrationPoint(
            name="s",
            age_constraint=SkewTConstraint(
                location=100.0, scale=10.0, shape=0.5, df=10.0
            ),
            resolved_taxa=["C", "D"],
            mrca_leaf_pair=("C", "D"),
        )

        summary = ac.get_calibration_summary([cal_gamma, cal_skew])

        # Gamma/Skew 约束必须走 describe()，不再显示 unknown
        assert "unknown" not in summary
        assert "Gamma" in summary
        assert "SkewT" in summary


# --------------------------------------------------------------------------- #
# 13) mcmctree_method：RootAge 单位已修正为 Ga（+ DEFAULT_ROOT_AGE_GA 常量）
# --------------------------------------------------------------------------- #
class TestMcmctreeRootAgeGa:
    def test_default_root_age_ga_constant(self):
        from phylodater.adapters.mcmctree_method import DEFAULT_ROOT_AGE_GA

        assert DEFAULT_ROOT_AGE_GA == 10.0
        assert isinstance(DEFAULT_ROOT_AGE_GA, float)

    def test_root_age_converted_to_ga_in_control_file(self, tmp_path):
        from phylodater.adapters.mcmctree_method import MCMCTreeMethod

        config = SimpleNamespace(
            rate_alpha=1.0,
            model="auto",
            ndata=1,
            clock=1,
            clean_data=0,
            burnin=2000,
            sampfreq=100,
            nsample=20000,
            root_age=None,
            paml_version="4.10.8",
            time_unit="Ma",
        )
        common = SimpleNamespace(seed=None)
        method = MCMCTreeMethod(config, output_dir=tmp_path, common_config=common)

        # 根节点校准：fixed_age = 3500 Ma
        root_cal = CalibrationPoint(
            name="ROOT_TEST",
            age_constraint=FixedAgeConstraint(fixed_age=3500.0),
            resolved_taxa=["A", "B"],
            mrca_leaf_pair=("A", "B"),
            is_root_node=True,
        )
        method._calibrations = [root_cal]

        ctl = tmp_path / "mcmctree.ctl"
        # 仅生成控制文件内容，不执行任何外部二进制
        method._generate_control_file(ctl, usedata=3)

        content = ctl.read_text(encoding="utf-8")
        rootage_lines = [
            ln for ln in content.splitlines() if ln.strip().startswith("RootAge")
        ]
        assert rootage_lines, "RootAge line not found in generated control file"

        # 3500 Ma -> 3.5 Ga（修复点：此前未 ÷1000，会得到 '<3500'）
        assert "RootAge = '<3.5'" in rootage_lines[0], rootage_lines[0]
