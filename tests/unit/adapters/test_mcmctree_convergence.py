"""
MCMCTree 收敛诊断单元测试

两部分：

1. 纯数学函数（wave-1 #8 修复验证）：标准 Gelman-Rubin PSRF 与 Geyer (1992)
   initial positive sequence ESS 的数值行为——收敛链 PSRF ≈ 1、白噪声链
   ESS ≈ n、自相关链 ESS 明显降低、均值不同的链 PSRF ≫ 1。
2. **真实 mcmc.txt 格式**的回归测试（审阅报告 A-1 / A-2 / B-4）：
   fixture 的表头按上游 ``collectx()`` 逐字生成
   （``mcmctree.c`` 4.10.8:2045-2110、4.8:1724-1741）::

       Gen\\tt_n5\\tt_n6\\tt_n7\\tmu\\tsigma2\\tlnL

   第一列 ``Gen`` 是非数值表头 + 各链完全相同的迭代号，正是旧实现
   ``np.loadtxt(mcmc_file)`` 抛 ValueError（A-1）与把正常运行误判为
   「停滞链」（A-2）的根源。这些测试断言的是 ``_check_convergence()`` 的
   文件加载 + 判定分支，而不是只测两个纯函数。
"""

import math
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest

from phylodater.adapters.mcmctree_method import (
    DEFAULT_ESS_THRESHOLD,
    DEFAULT_PSRF_THRESHOLD,
    MCMCTreeMethod,
    _diagnose_mcmc_chains,
    _gelman_rubin_psrf,
    _geyer_ess,
    _read_mcmc_chain,
    _split_iteration_columns,
    _summarize_chain_ages,
)
from phylodater.infrastructure.configuration import CommonConfig, MCMCTreeConfig

# --------------------------------------------------------------------------- #
# 真实格式 fixture 构造
# --------------------------------------------------------------------------- #
# 4 叶树（nnode=7）单基因（g=1）clock>1 的 mcmctree 表头，逐字取自上游：
#   fprintf "Gen" + "\tt_n%d" (i=s..nnode-1) + "\tmu" + "\tsigma2" + "\tlnL"
MCMC_HEADER = "Gen\tt_n5\tt_n6\tt_n7\tmu\tsigma2\tlnL"
MCMC_PARAMETER_COLUMNS = ["t_n5", "t_n6", "t_n7", "mu", "sigma2", "lnL"]


def _normal_chains(m: int, n: int, seed: int):
    rng = np.random.default_rng(seed)
    return [rng.normal(0.0, 1.0, size=n) for _ in range(m)]


def _write_mcmc_txt(
    path: Path,
    columns,
    n_samples: int,
    sampfreq: int = 1000,
    header: str = MCMC_HEADER,
) -> None:
    """按 MCMCTree 的真实排版写出 mcmc.txt（表头 + 制表符分隔的采样行）。

    Args:
        columns: 与表头中非 ``Gen`` 列一一对应的数组列表（长度 n_samples）。
    """
    header_columns = header.split("\t")
    assert len(header_columns) == len(columns) + 1, "列数与表头不匹配"
    lines = [header]
    for i in range(n_samples):
        fields = [str((i + 1) * sampfreq)]  # 迭代号列（各链逐元素相同）
        for values in columns:
            fields.append(f"{float(values[i]):.7f}")
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _method_with_chain_files(tmp_path: Path, chains_columns, num_runs: int, **cfg_kw):
    """创建适配器并把每条链写成 work_dir/run{rid}/mcmc.txt。"""
    config = MCMCTreeConfig(num_runs=num_runs, **cfg_kw)
    method = MCMCTreeMethod(
        config, output_dir=tmp_path, common_config=CommonConfig(nthreads=1)
    )
    for run_id, columns in enumerate(chains_columns, 1):
        run_dir = method.work_dir / f"run{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)
        n_samples = len(columns[0])
        _write_mcmc_txt(run_dir / "mcmc.txt", columns, n_samples)
    return method


def _iid_columns(n: int, seed: int, scale: float = 0.02, offset: float = 0.0):
    """生成一组（近似独立、无自相关）的同分布参数列。"""
    rng = np.random.default_rng(seed)
    return [
        offset + 0.30 + scale * rng.normal(size=n),  # t_n5
        offset + 0.22 + scale * rng.normal(size=n),  # t_n6
        offset + 0.15 + scale * rng.normal(size=n),  # t_n7
        0.5 + 0.05 * rng.normal(size=n),  # mu
        0.3 + 0.05 * rng.normal(size=n),  # sigma2
        -100.0 + 1.0 * rng.normal(size=n),  # lnL
    ]


def _ar1_columns(n: int, seed: int, rho: float):
    """强自相关（AR1）参数列：同分布但 ESS 远低于样本数。"""
    rng = np.random.default_rng(seed)

    def ar1(mean, sd):
        noise = rng.normal(0.0, sd * math.sqrt(1 - rho**2), size=n)
        out = np.zeros(n)
        out[0] = mean + rng.normal(0.0, sd)
        for i in range(1, n):
            out[i] = mean + rho * (out[i - 1] - mean) + noise[i]
        return out

    return [
        ar1(0.30, 0.02),
        ar1(0.22, 0.02),
        ar1(0.15, 0.02),
        ar1(0.5, 0.05),
        ar1(0.3, 0.05),
        ar1(-100.0, 1.0),
    ]


def _sawtooth_columns(n: int, shift: float):
    """确定性参数列：组内方差已知、链间均值偏移 = ``shift``。

    交替 ±1 序列的组内方差 ≈ 1、链间方差 = shift²/2，因此
    PSRF ≈ sqrt(1 + shift²/2)（与随机噪声不同，可精确复现）：
    shift=0.0 → PSRF ≈ 1.0；shift=0.7 → PSRF ≈ 1.12（旧实现 1.2 阈值下判「已收敛」）。
    该序列自相关为 -1，`_geyer_ess` 视作近乎独立 → ESS 很大，从而把「ESS 不足」
    与「PSRF 过大」两条判定分离开。
    """
    base = np.tile([-1.0, 1.0], n // 2 + 1)[:n]
    return [
        base + shift,  # t_n5：唯一有链间偏移的列
        base * 0.5 + 0.3,  # t_n6
        base * 0.2 + 0.15,  # t_n7
        base * 0.05 + 0.5,  # mu
        base * 0.05 + 0.3,  # sigma2
        base * 0.02 - 100.0,  # lnL
    ]


# --------------------------------------------------------------------------- #
# 1) 纯数学函数（原有覆盖，保持不变）
# --------------------------------------------------------------------------- #
def test_psrf_converged_chains_close_to_one():
    """同分布多链应给出 PSRF ≈ 1。"""
    chains = _normal_chains(4, 2000, seed=1)
    psrf, zero_var = _gelman_rubin_psrf(chains)
    assert not zero_var
    assert psrf == pytest.approx(1.0, abs=0.1)


def test_psrf_detects_nonconvergence():
    """均值显著不同的链应给出很大的 PSRF。"""
    rng = np.random.default_rng(4)
    chains = [rng.normal(0.0, 1.0, 2000), rng.normal(10.0, 1.0, 2000)]
    psrf, zero_var = _gelman_rubin_psrf(chains)
    assert not zero_var
    assert psrf > 2.0


def test_psrf_zero_variance_chain():
    """恒定（零方差）链标记为停滞（is_zero_variance=True, PSRF=inf）。"""
    chains = [np.full(500, 3.0), np.full(500, 3.0)]
    psrf, zero_var = _gelman_rubin_psrf(chains)
    assert zero_var
    assert psrf == float("inf")


def test_ess_white_noise_approx_n():
    """白噪声链的有效样本量应接近样本数 n。"""
    rng = np.random.default_rng(2)
    chain = rng.normal(0.0, 1.0, size=2000)
    ess = _geyer_ess(chain)
    n = len(chain)
    assert ess > 0.5 * n
    assert ess < 1.5 * n


def test_ess_autocorrelation_reduces_effective_size():
    """强自相关的 AR(1) 链有效样本量应明显低于白噪声。"""
    n = 2000
    rng = np.random.default_rng(3)
    noise = rng.normal(0.0, 1.0, size=n)
    rho = 0.5
    ar1 = np.zeros(n)
    for i in range(1, n):
        ar1[i] = rho * ar1[i - 1] + noise[i]

    ess_ar1 = _geyer_ess(ar1)
    ess_white = _geyer_ess(noise)

    # 自相关降低有效样本量
    assert ess_ar1 > 0
    assert ess_ar1 < 0.8 * ess_white
    # 理论 ESS ≈ n * (1 - rho) / (1 + rho) = n / 3 ≈ 666
    assert ess_ar1 > 100


def test_ess_constant_chain_is_zero():
    """恒定链方差为 0 → ESS = 0。"""
    chain = np.full(1000, 5.0)
    assert _geyer_ess(chain) == 0.0


# --------------------------------------------------------------------------- #
# 2) A-1：真实 mcmc.txt 必须能被加载（旧实现第一条链就抛 ValueError）
# --------------------------------------------------------------------------- #
def test_real_mcmc_txt_rejects_bare_loadtxt(tmp_path):
    """复现 A-1 的病因：np.loadtxt 无法直接读 MCMCTree 的表头行。"""
    columns = _iid_columns(50, seed=11)
    path = tmp_path / "mcmc.txt"
    _write_mcmc_txt(path, columns, 50)

    with pytest.raises(ValueError):
        np.loadtxt(path)  # 旧实现的做法


def test_read_mcmc_chain_drops_iteration_column_by_header_name(tmp_path):
    """_read_mcmc_chain 跳过表头并按列名剔除迭代号列（A-1 + A-2）。"""
    columns = _iid_columns(80, seed=12)
    path = tmp_path / "mcmc.txt"
    _write_mcmc_txt(path, columns, 80)

    samples, names = _read_mcmc_chain(path)

    assert names == MCMC_PARAMETER_COLUMNS
    assert samples.shape == (80, 6)
    # 迭代号列已被剔除：剩下的每一列都有非零方差
    assert all(
        float(np.var(samples[:, i], ddof=1)) > 0 for i in range(samples.shape[1])
    )


def test_split_iteration_columns_is_case_insensitive_and_name_based():
    """迭代号列按列名（不分大小写、不看位置）识别。"""
    keep, dropped = _split_iteration_columns(["Gen", "t_n5", "mu", "lnL"])
    assert keep == [1, 2, 3]
    assert dropped == ["Gen"]

    keep, dropped = _split_iteration_columns(["t_n5", "ITER", "mu"])
    assert keep == [0, 2]
    assert dropped == ["ITER"]


def test_read_mcmc_chain_without_header_still_parses(tmp_path):
    """无表头（被截断/异常输出）时按位置命名，不丢样本行。"""
    path = tmp_path / "mcmc.txt"
    rows = [
        "\t".join([str(1000 * (i + 1)), f"{0.3 + 0.01 * i:.7f}", f"{-100 + i:.3f}"])
        for i in range(6)
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    samples, names = _read_mcmc_chain(path)
    assert names == ["col1", "col2", "col3"]
    assert samples.shape == (6, 3)


def test_read_mcmc_chain_rejects_header_only_file(tmp_path):
    """只有表头、没有任何采样的 mcmc.txt 必须报错（而不是静默判收敛）。"""
    path = tmp_path / "mcmc.txt"
    path.write_text(MCMC_HEADER + "\n", encoding="utf-8")

    with pytest.raises(ValueError):
        _read_mcmc_chain(path)


def test_read_mcmc_chain_rejects_header_column_mismatch(tmp_path):
    """表头列数与采样列数不一致时报错。"""
    path = tmp_path / "mcmc.txt"
    path.write_text(MCMC_HEADER + "\n1000\t0.3000000\t0.2000000\n", encoding="utf-8")
    with pytest.raises(ValueError):
        _read_mcmc_chain(path)


# --------------------------------------------------------------------------- #
# 3) A-2：迭代号列不得被当成「停滞链」
# --------------------------------------------------------------------------- #
def test_gen_column_is_excluded_from_diagnostics():
    """A-2：迭代号列不得参与 PSRF/ESS——它的 ESS 极低，会拖垮整次运行的判定。"""
    n = 300
    ramp = np.arange(1, n + 1) * 1000.0
    chains = []
    for run in range(2):
        columns = _iid_columns(n, seed=20 + run)
        chains.append(np.column_stack([ramp] + columns))
    names = ["Gen"] + MCMC_PARAMETER_COLUMNS

    # 病因演示：单看迭代号列，它的自相关时间极长 → ESS 远低于随机列
    assert _geyer_ess(ramp) < DEFAULT_ESS_THRESHOLD
    assert _geyer_ess(chains[0][:, 1]) > DEFAULT_ESS_THRESHOLD

    result = _diagnose_mcmc_chains(chains, names)
    assert result["index_columns"] == ["Gen"]
    assert result["zero_variance_columns"] == []
    # 剔除后 min_ess 只由随机参数列决定，正常运行不会被迭代号列判死
    assert result["min_ess"] > DEFAULT_ESS_THRESHOLD
    assert result["max_psrf"] == pytest.approx(1.0, abs=0.05)


def test_constant_stochastic_column_is_not_mistaken_for_index():
    """零方差且非单调的列是真停滞，必须留在 zero_variance_columns 里。"""
    n = 200
    columns_a = _iid_columns(n, seed=23)
    columns_b = _iid_columns(n, seed=24)
    for columns in (columns_a, columns_b):
        columns[3] = np.full(n, 0.5)
    chains = [np.column_stack(columns) for columns in (columns_a, columns_b)]

    result = _diagnose_mcmc_chains(chains, MCMC_PARAMETER_COLUMNS)
    assert result["index_columns"] == []
    assert result["zero_variance_columns"] == ["mu"]


def test_check_convergence_real_format_reports_converged(tmp_path):
    """真实格式、三条同分布链 → 判为已收敛（A-1 修复前恒为 False）。"""
    n = 400
    chains = [
        _iid_columns(n, seed=31),
        _iid_columns(n, seed=32),
        _iid_columns(n, seed=33),
    ]
    method = _method_with_chain_files(tmp_path, chains, num_runs=3)
    recorder = Mock()
    method.logger = recorder

    assert method._check_convergence() is True

    diagnostics = method._convergence_diagnostics
    assert diagnostics["n_chains"] == 3
    assert diagnostics["n_samples_per_chain"] == n
    assert diagnostics["zero_variance_columns"] == []
    assert diagnostics["max_psrf"] <= diagnostics["psrf_threshold"]
    assert diagnostics["min_ess"] >= diagnostics["ess_threshold"]
    # 绝不能再出现「加载失败」或「停滞链」的误诊文案
    logged = " ".join(str(call) for call in recorder.method_calls)
    assert "Failed to load mcmc.txt" not in logged
    assert "stuck chains" not in logged


def test_check_convergence_logs_no_warning_for_real_format(tmp_path):
    """真实格式下不应触发 numpy/加载类告警。"""
    n = 300
    chains = [_iid_columns(n, seed=41), _iid_columns(n, seed=42)]
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)
    recorder = Mock()
    method.logger = recorder

    assert method._check_convergence() is True
    warnings = " ".join(str(c.args) for c in recorder.warning.call_args_list)
    assert "mcmc.txt" not in warnings


# --------------------------------------------------------------------------- #
# 4) B-4：阈值必须把明显未收敛判为未收敛
# --------------------------------------------------------------------------- #
def test_thresholds_are_tightened_defaults():
    """默认阈值不再是 1.2 / 2.0，且 ESS 也参与布尔判定。"""
    assert DEFAULT_PSRF_THRESHOLD <= 1.1
    assert DEFAULT_ESS_THRESHOLD >= 100


def test_psrf_between_1_05_and_1_2_is_not_converged(tmp_path):
    """旧实现「PSRF < 1.2 且 ESS 达标 → True」这一档必须改判未收敛。"""
    n = 400
    # shift=0.7 → PSRF ≈ sqrt(1 + 0.7²/2) ≈ 1.12：> 新阈值 1.05，< 旧阈值 1.2
    chains = [_sawtooth_columns(n, shift=0.0), _sawtooth_columns(n, shift=0.7)]
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)
    recorder = Mock()
    method.logger = recorder

    diagnostics_before = method._check_convergence()
    diagnostics = method._convergence_diagnostics
    assert diagnostics_before is False
    assert 1.05 < diagnostics["max_psrf"] < 1.2
    assert diagnostics["min_ess"] >= diagnostics["ess_threshold"]
    assert diagnostics["verdict"] == "not_converged"
    errors = " ".join(
        " ".join(str(a) for a in c.args) for c in recorder.error.call_args_list
    )
    assert "NOT converged" in errors
    # 日志文案必须与返回值一致：旧实现写「may not have fully converged」却 return True
    assert "may not have fully converged" not in errors


def test_low_ess_with_good_psrf_is_not_converged(tmp_path):
    """PSRF 达标但 ESS 不足时不得判收敛（旧实现 return True）。"""
    n = 500
    chains = [_ar1_columns(n, seed=61, rho=0.985), _ar1_columns(n, seed=62, rho=0.985)]
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)
    recorder = Mock()
    method.logger = recorder

    assert method._check_convergence() is False
    diagnostics = method._convergence_diagnostics
    assert diagnostics["min_ess"] < diagnostics["ess_threshold"]
    assert diagnostics["verdict"] == "not_converged"
    assert any(
        "ESS" in " ".join(str(a) for a in c.args) for c in recorder.error.call_args_list
    )


def test_ess_only_gate_is_enforced_via_config(tmp_path):
    """独立同分布链（PSRF≈1）在 ESS 阈值更高时必须判未收敛——两阈值是「与」关系。"""
    n = 400
    chains = [_sawtooth_columns(n, shift=0.0), _sawtooth_columns(n, shift=0.0)]
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)

    assert method._check_convergence() is True
    # 交替序列近乎独立，ESS 被算得极大；把阈值抬到它之上即可单独检验 ESS 分支
    method.config.ess_threshold = 10**15
    assert method._check_convergence() is False
    diagnostics = method._convergence_diagnostics
    assert diagnostics["max_psrf"] <= diagnostics["psrf_threshold"]
    assert diagnostics["min_ess"] < diagnostics["ess_threshold"]
    assert any("ESS" in reason for reason in diagnostics["reasons"])


def test_thresholds_can_be_overridden_by_config(tmp_path):
    """阈值可被配置项覆盖（不再是只藏在代码里的魔法数字）。"""
    n = 400
    chains = [_sawtooth_columns(n, shift=0.7), _sawtooth_columns(n, shift=0.0)]
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)

    assert method._check_convergence() is False
    default_psrf = method._convergence_diagnostics["psrf_threshold"]
    assert default_psrf == pytest.approx(DEFAULT_PSRF_THRESHOLD)

    method.config.psrf_threshold = 1.5
    assert method._check_convergence() is True
    assert method._convergence_diagnostics["psrf_threshold"] == 1.5


def test_stuck_stochastic_parameter_reports_column_name(tmp_path):
    """真正停滞（零方差）的随机参数列必须判未收敛，并点名其列名。"""
    n = 300
    chains = [_iid_columns(n, seed=81), _iid_columns(n, seed=82)]
    for columns in chains:
        columns[3] = np.full(n, 0.5)  # mu 完全不动
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)
    recorder = Mock()
    method.logger = recorder

    assert method._check_convergence() is False
    assert method._convergence_diagnostics["zero_variance_columns"] == ["mu"]
    assert any(
        "stuck chains" in " ".join(str(a) for a in c.args)
        for c in recorder.error.call_args_list
    )


def test_nan_samples_are_not_converged(tmp_path):
    """链中含 NaN → 判未收敛（不得静默跳过该列）。"""
    n = 200
    chains = [_iid_columns(n, seed=91), _iid_columns(n, seed=92)]
    chains[1][1][50] = float("nan")
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)

    assert method._check_convergence() is False


def test_missing_mcmc_txt_is_not_converged(tmp_path):
    """num_runs>=2 但缺 mcmc.txt：无法核验 → fail-closed 判未收敛。"""
    config = MCMCTreeConfig(num_runs=2)
    method = MCMCTreeMethod(config, output_dir=tmp_path, common_config=CommonConfig())
    recorder = Mock()
    method.logger = recorder

    assert method._check_convergence() is False
    assert any(
        "mcmc.txt not found" in str(c.args) for c in recorder.warning.call_args_list
    )


def test_column_layout_mismatch_between_runs_is_not_diagnosable(tmp_path):
    """各链列不可对齐时返回 None（无法诊断），而不是 True/False。"""
    n = 200
    chains = [_iid_columns(n, seed=101), _iid_columns(n, seed=102)]
    method = _method_with_chain_files(tmp_path, chains, num_runs=2)
    extra = chains[1] + [np.zeros(n)]
    second_dir = method.work_dir / "run2"
    _write_mcmc_txt(
        second_dir / "mcmc.txt",
        extra,
        n,
        header=MCMC_HEADER + "\tmisc",
    )
    recorder = Mock()
    method.logger = recorder

    assert method._check_convergence() is None
    assert any(
        "column layout differs" in str(c.args) for c in recorder.error.call_args_list
    )


def test_single_run_returns_none(tmp_path):
    """单链运行不判收敛（返回 None）。"""
    config = MCMCTreeConfig(num_runs=1)
    method = MCMCTreeMethod(config, output_dir=tmp_path, common_config=CommonConfig())
    assert method._check_convergence() is None


# --------------------------------------------------------------------------- #
# 5) B-3：多链汇总的区间口径标注
# --------------------------------------------------------------------------- #
def test_summarize_chain_ages_single_chain_is_hpd():
    """单链区间就是 MCMCTree 自己的 95% HPD。"""
    from phylodater.models import CIType

    summary = _summarize_chain_ages(
        [{"mean": 0.2, "median": 0.19, "lower": 0.1, "upper": 0.3}]
    )
    assert summary["ci_type"] is CIType.HPD95
    assert summary["is_envelope"] is False
    assert summary["mean"] == pytest.approx(0.2)
    assert summary["median"] == pytest.approx(0.19)


def test_summarize_chain_ages_multi_chain_is_envelope_not_hpd():
    """多链区间是外包络 → 标 RANGE，并且 mean/median 来自同一批链。"""
    from phylodater.models import CIType

    summary = _summarize_chain_ages(
        [
            {"mean": 0.2, "median": 0.2, "lower": 0.1, "upper": 0.3},
            {"mean": 0.4, "median": 0.4, "lower": 0.35, "upper": 0.5},
        ]
    )
    assert summary["ci_type"] is CIType.RANGE
    assert summary["is_envelope"] is True
    assert summary["ci_lower"] == pytest.approx(0.1)
    assert summary["ci_upper"] == pytest.approx(0.5)
    # mean 与 median 都对两条链取，不再出现「mean 全链、median 只 run1」
    assert summary["mean"] == pytest.approx(0.3)
    assert summary["median"] == pytest.approx(0.3)


def test_summarize_chain_ages_detects_disjoint_chain_intervals():
    """两条链的 95% 区间完全不相交时必须暴露。"""
    summary = _summarize_chain_ages(
        [
            {"mean": 0.2, "median": 0.2, "lower": 0.1, "upper": 0.25},
            {"mean": 0.6, "median": 0.6, "lower": 0.45, "upper": 0.7},
        ]
    )
    assert summary["disjoint_chain_pairs"] == [(1, 2)]
