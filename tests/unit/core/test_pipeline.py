"""Per-method sub-seed forwarding (review items P0-1 / E-3).

Two distinct defects are pinned here.

P0-1 (historical): ``SeedSequence(base).spawn(n)`` children all report
``.entropy == base``, so reading ``.entropy`` collapsed every parallel method
onto one shared seed. ``_derive_method_seeds`` must draw each child's own value
so the seeds are distinct yet reproducible.

E-3 (fixed now): the derived seeds were only fed to the *process-level*
``random.seed()`` / ``np.random.seed()`` of each worker thread, while every
adapter reads the shared ``ToolConfig.common.seed``. All seven engines therefore
received the same seed. The sub-seed must travel **inside the configuration
object** handed to the adapter, and the worker must stop mutating global RNG
state.
"""

import random
import warnings as _warnings
from pathlib import Path

import numpy as np
import pytest

from phylodater.core.exceptions import SemanticDegradationWarning
from phylodater.core.method_interface import DatingMethod, DatingMethodRegistry
from phylodater.core.pipeline import (
    DatingPipeline,
    ParallelPipeline,
    PipelineConfig,
    _derive_method_seeds,
)
from phylodater.infrastructure.configuration import ToolConfig, derive_method_seed
from phylodater.models import DatingResult

ALL_METHODS = [
    "mcmctree",
    "lsd2",
    "treepl",
    "pathd8",
    "r8s",
    "wlogdate",
    "mdcat",
]

PROBE_NAMES = ["probe_a", "probe_b", "probe_c"]


class TestSubSeedDerivation:
    """``_derive_method_seeds``：互异、可复现、与本批方法数无关。"""

    def test_every_method_is_assigned(self):
        seeds = _derive_method_seeds(42, ALL_METHODS)
        assert set(seeds) == set(ALL_METHODS)

    def test_sub_seeds_are_pairwise_distinct(self):
        """The core P0-1 assertion: no two of the seven methods share a seed."""
        seeds = _derive_method_seeds(42, ALL_METHODS)
        values = list(seeds.values())
        assert len(set(values)) == len(
            values
        ), "P0-1 regression: parallel methods collapsed onto a shared seed"

    def test_sub_seeds_are_reproducible_for_fixed_base(self):
        a = _derive_method_seeds(42, ALL_METHODS)
        b = _derive_method_seeds(42, ALL_METHODS)
        assert a == b

    def test_different_base_seed_changes_streams(self):
        a = _derive_method_seeds(42, ["mcmctree", "lsd2"])
        b = _derive_method_seeds(1234, ["mcmctree", "lsd2"])
        assert a != b

    def test_single_and_empty_method_lists(self):
        one = _derive_method_seeds(7, ["mcmctree"])
        assert list(one) == ["mcmctree"]
        assert isinstance(one["mcmctree"], int)
        assert _derive_method_seeds(7, []) == {}

    def test_derivation_agrees_with_configuration_helper(self):
        """批量派生与逐名派生必须是同一张表（否则串行/并行会分叉）。"""
        seeds = _derive_method_seeds(99, ALL_METHODS)
        assert seeds == {name: derive_method_seed(99, name) for name in ALL_METHODS}

    def test_sub_seed_does_not_depend_on_batch_size(self):
        """E-3：串行 ``-j 1`` 与七方法并行下，同一方法拿到同一个子种子。"""
        solo = _derive_method_seeds(42, ["mcmctree"])["mcmctree"]
        full = _derive_method_seeds(42, ALL_METHODS)["mcmctree"]
        assert solo == full

    def test_unknown_base_seed_still_yields_distinct_seeds(self):
        """``base_seed=None`` 时也必须互异（同一次运行内可复盘）。"""
        seeds = _derive_method_seeds(None, ALL_METHODS)
        assert len(set(seeds.values())) == len(ALL_METHODS)


class SeedRecordingAdapter(DatingMethod):
    """记录流水线 handed 给本适配器的种子与随机流。"""

    #: 每个实例登记自己，供"执行之后回看"的断言使用
    instances: list = []

    def __init__(
        self, config, output_dir, software_paths=None, common_config=None
    ) -> None:
        super().__init__(
            config, output_dir, software_paths, common_config=common_config
        )
        self.recorded_seed = self.common_config.seed
        self.rng_at_creation = self.rng
        type(self).instances.append(self)

    @property
    def method_name(self) -> str:
        return "seed_probe"

    def validate_environment(self) -> bool:
        return True

    def prepare_inputs(self, tree, calibrations, alignment_path=None):
        return {}

    def execute(self) -> bool:
        return True

    def parse_results(self) -> DatingResult:
        return DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=";",
        )


@pytest.fixture
def probe_registry():
    """把三个探针名字注册到注册表，测完原样还回。"""
    saved = dict(DatingMethodRegistry._registry)
    SeedRecordingAdapter.instances = []
    for name in PROBE_NAMES:
        DatingMethodRegistry.register(name, SeedRecordingAdapter)
    try:
        yield
    finally:
        DatingMethodRegistry._registry.clear()
        DatingMethodRegistry._registry.update(saved)


def _pipeline(tmp_path, cls=DatingPipeline, seed=42, methods=("probe_a",)):
    config = PipelineConfig(
        methods=list(methods), output_dir=tmp_path, seed=seed, enable_checkpoint=False
    )
    return cls(config, tool_config=ToolConfig())


class TestMethodSeedForwarding:
    """E-3 核心：子种子随配置到达适配器，而不是只喂全局 RNG。"""

    def test_each_method_gets_its_own_sub_seed(self, tmp_path, probe_registry):
        pipeline = _pipeline(tmp_path, seed=42)

        adapters = {name: pipeline._create_adapter(name) for name in PROBE_NAMES}

        seen = [adapter.recorded_seed for adapter in adapters.values()]
        assert len(set(seen)) == len(PROBE_NAMES), "适配器仍拿到同一个种子"
        for name, adapter in adapters.items():
            assert adapter.recorded_seed == derive_method_seed(42, name)

    def test_shared_tool_config_is_not_mutated(self, tmp_path, probe_registry):
        pipeline = _pipeline(tmp_path, seed=42)

        pipeline._create_adapter("probe_a")
        pipeline._create_adapter("probe_b")

        assert pipeline._tool_config.common.seed == 42

    def test_adapter_gets_its_own_generator(self, tmp_path, probe_registry):
        pipeline = _pipeline(tmp_path, seed=42)

        first = pipeline._create_adapter("probe_a")
        second = pipeline._create_adapter("probe_b")

        assert isinstance(first.rng, np.random.Generator)
        assert first.rng is not second.rng
        assert first.method_seed != second.method_seed
        # 同一子种子的两条流必然给出同一串数（可复现）
        assert first.rng.random(3).tolist() == pipeline._create_adapter(
            "probe_a"
        ).rng.random(3).tolist()

    def test_same_base_seed_reproduces_across_pipelines(
        self, tmp_path, probe_registry
    ):
        a = _pipeline(tmp_path / "a", seed=1234)
        b = _pipeline(tmp_path / "b", seed=1234)

        assert [
            a._create_adapter(name).recorded_seed for name in PROBE_NAMES
        ] == [b._create_adapter(name).recorded_seed for name in PROBE_NAMES]

    def test_serial_execution_forwards_sub_seed(
        self, tmp_path, probe_registry, monkeypatch
    ):
        """串行路径同样必须走方法级种子（不是一条只在并行才生效的路径）。"""
        pipeline = _pipeline(tmp_path, seed=7, methods=PROBE_NAMES)
        monkeypatch.setattr(
            pipeline, "_plan_method_execution", lambda *a, **k: list(PROBE_NAMES)
        )

        pipeline._execute_methods(None, None, [])

        assert set(pipeline.results) == set(PROBE_NAMES)
        assert sorted(a.recorded_seed for a in SeedRecordingAdapter.instances) == sorted(
            derive_method_seed(7, name) for name in PROBE_NAMES
        )

    def test_parallel_execution_forwards_sub_seed(
        self, tmp_path, probe_registry, monkeypatch
    ):
        pipeline = _pipeline(
            tmp_path, cls=ParallelPipeline, seed=7, methods=PROBE_NAMES
        )
        monkeypatch.setattr(
            pipeline, "_plan_method_execution", lambda *a, **k: list(PROBE_NAMES)
        )

        pipeline._execute_methods(None, None, [])

        assert set(pipeline.results) == set(PROBE_NAMES)
        assert sorted(a.recorded_seed for a in SeedRecordingAdapter.instances) == sorted(
            derive_method_seed(7, name) for name in PROBE_NAMES
        )

    def test_parallel_worker_does_not_touch_global_rng(
        self, tmp_path, probe_registry, monkeypatch
    ):
        """E-3：worker 内不再给进程级 ``random`` / ``np.random`` 播种。"""
        pipeline = _pipeline(
            tmp_path, cls=ParallelPipeline, seed=7, methods=PROBE_NAMES
        )
        monkeypatch.setattr(
            pipeline, "_plan_method_execution", lambda *a, **k: list(PROBE_NAMES)
        )

        before_python = random.getstate()
        before_numpy = np.random.get_state()

        pipeline._execute_methods(None, None, [])

        assert set(pipeline.results) == set(
            PROBE_NAMES
        ), "方法执行失败（全局 RNG 若仍被播种，此处会捕获到异常）"
        assert random.getstate() == before_python, "全局 random 被 worker 改写了"
        after_numpy = np.random.get_state()
        assert after_numpy[0] == before_numpy[0] and np.array_equal(
            after_numpy[1], before_numpy[1]
        ), "全局 np.random 被 worker 改写了"


class TestRealEnginesGetOwnSeed:
    """四个"把种子写进引擎"的适配器必须各拿各的子种子（E-3 端到端证据）。

    ``mcmctree`` 写控制文件 ``seed =``、``treepl`` 写 ``seed =``、``wlogdate``
    传 ``-s``、``mdcat`` 传 ``--randSeed``——四处读的都是
    ``self.common_config.seed``，所以只要该字段逐方法不同，种子就真的转发了。
    """

    SEED_BEARING = ["mcmctree", "treepl", "wlogdate", "mdcat"]

    def test_registered_adapters_read_distinct_sub_seeds(self, tmp_path):
        import phylodater.adapters  # noqa: F401  触发注册

        pipeline = _pipeline(tmp_path, seed=42, methods=("mcmctree",))

        seeds = {
            name: pipeline._create_adapter(name).common_config.seed
            for name in self.SEED_BEARING
        }

        assert len(set(seeds.values())) == len(self.SEED_BEARING), (
            f"E-3 regression: 引擎仍共享同一个种子 {seeds}"
        )
        for name, seed in seeds.items():
            assert seed == derive_method_seed(42, name)

    def test_treepl_control_file_carries_method_seed(self, tmp_path):
        """子种子确实落到 treePL 控制文件文本里，而不只是内存对象上。"""
        import phylodater.adapters  # noqa: F401

        pipeline = _pipeline(tmp_path, seed=42, methods=("treepl",))
        adapter = pipeline._create_adapter("treepl")

        conf = tmp_path / "treepl_probe.conf"
        with open(conf, "w", encoding="utf-8") as handle:
            adapter._write_common_params(handle, Path("tree.nwk"))

        assert f"seed = {derive_method_seed(42, 'treepl')}" in conf.read_text(
            encoding="utf-8"
        )
        # 修复项 E-5：共享配置没有被就地改写
        assert pipeline._tool_config.common.seed == 42


class DegradationProbingAdapter(DatingMethod):
    """在 ``prepare_inputs`` 里触发三类**真实**语义降级，再交回带节点的结果。

    三类降级分别来自三个约束类、三个引擎，覆盖 E-8 要求（至少两类）：

    * ``UniformAgeConstraint`` → wLogDate：区间被压成中值点估计；
    * ``SoftLowerBoundConstraint`` → treePL：软下界被压成硬下界；
    * ``MaximumAgeConstraint`` → MD-Cat：仅上界被压成 ``max/2`` 点估计。
    """

    @property
    def method_name(self) -> str:
        return "degradation_probe"

    def validate_environment(self) -> bool:
        return True

    def prepare_inputs(self, tree, calibrations, alignment_path=None):
        from phylodater.models.constraints import (
            MaximumAgeConstraint,
            SoftLowerBoundConstraint,
            UniformAgeConstraint,
        )

        UniformAgeConstraint(min_age=2.0, max_age=4.0).to_software_format(
            "wlogdate", name="nodeA", taxa_str="A+B"
        )
        SoftLowerBoundConstraint(min_age=3.0).to_software_format(
            "treepl", name="nodeB"
        )
        MaximumAgeConstraint(max_age=5.0).to_software_format(
            "mdcat", name="nodeC", taxa_str="C+D"
        )
        return {"prepared": True}

    def execute(self) -> bool:
        return True

    def parse_results(self) -> DatingResult:
        from phylodater.models.results import NodeAgeEstimate

        return DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=";",
            node_ages={
                "nodeA": NodeAgeEstimate(mean_age=3.0),
                "nodeB": NodeAgeEstimate(mean_age=3.0),
                "nodeC": NodeAgeEstimate(mean_age=2.5),
            },
        )


@pytest.fixture
def degradation_registry():
    saved = dict(DatingMethodRegistry._registry)
    DatingMethodRegistry.register("degradation_probe", DegradationProbingAdapter)
    try:
        yield
    finally:
        DatingMethodRegistry._registry.clear()
        DatingMethodRegistry._registry.update(saved)


class TestDegradationWarningsArchived:
    """E-8：语义降级必须落进 ``DatingResult.warnings``；E-9：降级前形态可查。"""

    def _run(self, tmp_path):
        from phylodater.core.pipeline import run_adapter_and_capture_degradations

        adapter = DegradationProbingAdapter(ToolConfig(), tmp_path)
        return run_adapter_and_capture_degradations(adapter, None, [], None)

    def test_all_three_degradations_reach_result_warnings(self, tmp_path):
        result = self._run(tmp_path)
        texts = "\n".join(str(w) for w in result.warnings)
        assert len(result.warnings) == 3, result.warnings
        assert "wLogDate does not support range constraints" in texts
        assert "treepl does not support soft lower bounds" in texts
        assert "only supports point calibrations and does not " in texts

    def test_structured_records_carry_pre_degradation_shape(self, tmp_path):
        """E-9：结构化记录给出"降级前长什么样 / 实际交出去的是什么"。"""
        result = self._run(tmp_path)
        records = result.metadata["constraint_degradations"]
        assert len(records) == 3
        by_target = {r["target"]: r for r in records}
        assert set(by_target) == {"wlogdate", "treepl", "mdcat"}
        wlogdate = by_target["wlogdate"]
        assert wlogdate["constraint"] == "UniformAgeConstraint"
        # 规范化形态：逐字段（含带默认值的字段）打印，可机器比对
        assert wlogdate["original"] == (
            "UniformAgeConstraint(min_age=2.0, max_age=4.0, "
            "tail_lower=0.0001, tail_upper=0.0001)"
        )
        # 人类可读形态：与约束自己的 describe() 一致
        assert wlogdate["original_description"]
        assert "2.0-4.0 Ma" in wlogdate["original_description"]
        assert "3.0" in wlogdate["degraded_to"]
        assert wlogdate["node"] == "nodeA"
        assert wlogdate["message"]

    def test_affected_node_estimates_expose_degraded_from(self, tmp_path):
        result = self._run(tmp_path)
        degraded = result.node_ages["nodeA"]
        assert degraded.degraded_from is not None
        assert degraded.degraded_from["target"] == "wlogdate"
        assert (
            degraded.degraded_from["original"]
            == "UniformAgeConstraint(min_age=2.0, max_age=4.0, "
            "tail_lower=0.0001, tail_upper=0.0001)"
        )
        # 降级不可逆：年龄本身仍然是"降级后"的那个值
        assert degraded.mean_age == 3.0
        # 序列化里带得出去（审计时能从磁盘查回事前形态）
        assert degraded.to_dict()["degraded_from"]["target"] == "wlogdate"

    def test_warnings_deduplicated_and_adapter_entries_kept_first(self, tmp_path):
        """适配器自己写入的条目（r8s 的 C-41 形态）排在前面且不重复。"""
        from phylodater.core.pipeline import run_adapter_and_capture_degradations

        class PreSeededAdapter(DegradationProbingAdapter):
            def parse_results(self) -> DatingResult:
                result = super().parse_results()
                result.warnings.append("wLogDate does not support range constraints. "
                                       "Using midpoint 3.0 from range [2.0, 4.0]")
                return result

        adapter = PreSeededAdapter(ToolConfig(), tmp_path)
        result = run_adapter_and_capture_degradations(adapter, None, [], None)
        texts = [str(w) for w in result.warnings]
        assert len(texts) == len(set(texts)), "同一条降级被记了两遍"
        assert texts[0].startswith("wLogDate does not support range constraints")

    def test_non_degradation_warnings_are_not_swallowed(self, tmp_path):
        """窗口里收到的其它类别警告必须重发，不能因为落档而降级成"看不见"。"""
        import warnings as _warnings

        from phylodater.core.pipeline import run_adapter_and_capture_degradations

        class NoisyAdapter(DegradationProbingAdapter):
            def prepare_inputs(self, tree, calibrations, alignment_path=None):
                _warnings.warn("noise from elsewhere", UserWarning, stacklevel=2)
                return super().prepare_inputs(tree, calibrations, alignment_path)

        adapter = NoisyAdapter(ToolConfig(), tmp_path)
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            result = run_adapter_and_capture_degradations(adapter, None, [], None)
        assert not any("noise from elsewhere" in str(w) for w in result.warnings)
        assert any("noise from elsewhere" == str(w.message) for w in caught)

    def test_without_unified_exit_constraint_degradations_are_lost(
        self, tmp_path
    ):
        """E-8 的"修复前"锚点：绕过统一出口，约束级降级就只活在警告流里。

        这条测试钉住缺陷本身而不是巧合：三段调用照旧直跑，
        ``result.warnings`` 应当**为空**，而警告确实发出来了（3 条）。
        若哪天有人把 ``models/constraints.py`` 的降级改成直接写结果，
        这条会失败并提醒复核 E-8 的口径。
        """
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            adapter = DegradationProbingAdapter(ToolConfig(), tmp_path)
            adapter.prepare_inputs(None, [], None)
            adapter.execute()
            result = adapter.parse_results()

        assert len(caught) == 3
        assert all(
            issubclass(w.category, SemanticDegradationWarning) for w in caught
        )
        assert result.warnings == [], "约束级降级竟已自己进了结果，E-8 口径需重审"

    def test_serial_pipeline_path_archives_degradations(
        self, tmp_path, degradation_registry, monkeypatch
    ):
        """串行路径（``_execute_methods``）也走同一个统一出口。"""
        pipeline = _pipeline(tmp_path, methods=("degradation_probe",), seed=42)
        monkeypatch.setattr(
            pipeline, "_plan_method_execution", lambda *a, **k: ["degradation_probe"]
        )

        pipeline._execute_methods(None, None, [])

        result = pipeline.results["degradation_probe"]
        assert len(result.warnings) == 3
        assert len(result.metadata["constraint_degradations"]) == 3

    def test_parallel_pipeline_path_archives_degradations(
        self, tmp_path, degradation_registry, monkeypatch
    ):
        """并行路径同样落档 —— 口径不只在串行才成立。"""
        pipeline = _pipeline(
            tmp_path, cls=ParallelPipeline, methods=("degradation_probe",), seed=42
        )
        monkeypatch.setattr(
            pipeline, "_plan_method_execution", lambda *a, **k: ["degradation_probe"]
        )

        pipeline._execute_methods(None, None, [])

        result = pipeline.results["degradation_probe"]
        assert len(result.warnings) == 3
        assert len(result.metadata["constraint_degradations"]) == 3
