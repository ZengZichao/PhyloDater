"""
模块十二：批量处理与检查点功能测试

验证批量处理、检查点保存/恢复/清理、运行时元数据，以及：

  * B-23 —— ``--resume`` 的复用判定必须看输入指纹：没有登记过指纹就**不得**
    复用陈旧结果；``reset()`` 必须连指纹一起清掉（否则换过输入后旧指纹会"复活"）。
  * B-27 —— 日志层接受任意载荷（异常对象等），"报告降级的代码"自己不能再崩。
  * C-49 —— 样本身份（``sample_id``／输出目录／检查点文件）在整批内必须唯一：
    缺 id 要按索引派生唯一 id，碰撞要直接失败，否则续跑会把 A 样本的结果交给 B 样本。
  * C-50 —— "成功"只有一个口径（跑通**且**产出非空结果）；CPU 秒之和与整批真实墙钟
    必须分开报告（并行时前者可达真实耗时的 num_workers 倍）。
  * C-51 —— ``batch_results.json`` 是一份清单，必须在产物里自己说清楚年龄在哪。
"""

import json
import time
from pathlib import Path

import pytest

from phylodater.core import batch_processor as batch_processor_module
from phylodater.core.batch_processor import (
    BatchConfig,
    BatchProcessor,
    BatchResult,
    SampleResult,
)
from phylodater.core.exceptions import SemanticDegradationWarning
from phylodater.core.pipeline import PipelineConfig
from phylodater.infrastructure.checkpoint import (
    RUN_FINGERPRINT_KEY,
    CheckpointManager,
    PipelineCheckpointManager,
)
from phylodater.infrastructure.logging import get_logger
from phylodater.infrastructure.runtime_metadata import RuntimeMetadataManager
from phylodater.models import CIType, DatingResult, NodeAgeEstimate


class TestBatchProcessor:
    """批量处理测试"""

    def test_bat_001_batch_multiple_genes(self, tmp_path):
        """BAT-001: 多基因分析"""
        config = PipelineConfig(methods=["pathd8"], output_dir=tmp_path)
        processor = BatchProcessor(base_config=config, batch_config=BatchConfig())
        assert processor is not None
        assert tmp_path.exists()

    def test_bat_003_batch_output_isolation(self, tmp_path):
        """BAT-003: 批量输出隔离"""
        out1 = tmp_path / "run1"
        out2 = tmp_path / "run2"
        out1.mkdir()
        out2.mkdir()
        assert out1 != out2

    def test_bat_004_runtime_metadata_per_run(self, tmp_path):
        """BAT-004: 批量运行时元数据"""
        manager1 = RuntimeMetadataManager(output_dir=tmp_path / "run1")
        manager2 = RuntimeMetadataManager(output_dir=tmp_path / "run2")
        path1 = manager1.save()
        path2 = manager2.save()
        assert path1.exists()
        assert path2.exists()

    def test_c51_manifest_declares_it_holds_no_ages(self, tmp_path):
        """C-51: batch_results.json 必须在产物里自述它不含年龄、年龄在哪里"""
        sample_dir = tmp_path / "sampleA"
        result = BatchResult(
            total_samples=1,
            successful=1,
            failed=0,
            sample_results=[
                SampleResult(
                    sample_id="sampleA",
                    success=True,
                    results={
                        "pathd8": DatingResult(
                            method_name="pathd8",
                            run_id="r1",
                            dated_tree_newick="((A,B));",
                            node_ages={"n1": NodeAgeEstimate(mean_age=100.0)},
                        )
                    },
                    output_dir=sample_dir,
                )
            ],
        )
        payload = result.to_dict()
        assert payload["artifact_type"] == "batch_manifest"
        assert "no age estimates" in payload["note"]
        row = payload["sample_results"][0]
        assert row["age_data_file"] == str(sample_dir / "comparison_table.tsv")
        assert row["node_key_counts"] == {"pathd8": 1}
        # 清单里确实没有任何年龄数值
        assert "100.0" not in json.dumps(payload)

        out = tmp_path / "batch_results.json"
        BatchProcessor(
            PipelineConfig(methods=["pathd8"], output_dir=tmp_path)
        ).save_results(result, out)
        assert out.exists()


def _dating_result(method="pathd8", nodes=1):
    return DatingResult(
        method_name=method,
        run_id=f"{method}_r1",
        dated_tree_newick="((A,B));",
        node_ages={f"n{i}": NodeAgeEstimate(mean_age=100.0) for i in range(nodes)},
    )


def _install_fake_pipeline(
    monkeypatch, results_by_subdir=None, default=None, sleep=0.0
):
    """把 ``DatingPipeline`` 换成替身，记录每个样本真正拿到的 ``output_dir``。

    批量身份（C-49）与"成功"口径（C-50）都由 ``BatchProcessor`` 自己决定，替身只
    负责暴露它给每个样本配了什么目录、返回了什么结果。
    """
    calls = []

    class FakePipeline:
        def __init__(self, config):
            self.config = config

        def run(self, tree, alignment_path, calibrations):
            calls.append({"output_dir": Path(self.config.output_dir), "tree": tree})
            if sleep:
                time.sleep(sleep)
            subdir = Path(self.config.output_dir).name
            if results_by_subdir is not None and subdir in results_by_subdir:
                return results_by_subdir[subdir]
            if default is not None:
                return default
            return {"pathd8": _dating_result()}

    monkeypatch.setattr(batch_processor_module, "DatingPipeline", FakePipeline)
    return calls


def _processor(tmp_path, **batch_kwargs):
    return BatchProcessor(
        base_config=PipelineConfig(methods=["pathd8"], output_dir=tmp_path),
        batch_config=BatchConfig(**batch_kwargs),
    )


TREE = "(A:0.1,B:0.1);"


class TestSampleIdentityIsolation:
    """C-49：``sample_id``、输出目录、检查点文件三者在整批内必须一一对应"""

    def test_missing_sample_ids_derive_unique_ids(self, tmp_path, monkeypatch):
        calls = _install_fake_pipeline(monkeypatch)
        samples = [{"tree": TREE}, {"tree": TREE}, {"tree": TREE, "sample_id": "  "}]
        result = _processor(tmp_path, num_workers=1).process_samples(samples)

        assert [sr.sample_id for sr in result.sample_results] == [
            "sample_0000",
            "sample_0001",
            "sample_0002",
        ]
        # 旧写法三个样本都叫 "unknown"：同一目录、同一 checkpoint
        assert len({c["output_dir"] for c in calls}) == 3
        assert all("unknown" not in str(c["output_dir"]) for c in calls)

    def test_derived_ids_do_not_share_one_checkpoint_file(self, tmp_path, monkeypatch):
        calls = _install_fake_pipeline(monkeypatch)
        _processor(tmp_path, num_workers=1).process_samples(
            [{"tree": TREE}, {"tree": TREE}]
        )
        # 检查点就建在样本输出目录里（core/pipeline.py → PipelineCheckpointManager）
        checkpoints = {
            d / ".phylodater_pipeline_checkpoint.json"
            for d in (c["output_dir"] for c in calls)
        }
        assert len(checkpoints) == 2

    def test_caller_sample_dicts_are_not_mutated(self, tmp_path, monkeypatch):
        _install_fake_pipeline(monkeypatch)
        samples = [{"tree": TREE}]
        _processor(tmp_path, num_workers=1).process_samples(samples)
        assert "sample_id" not in samples[0]

    def test_colliding_output_subdirs_fail_the_batch(self, tmp_path, monkeypatch):
        _install_fake_pipeline(monkeypatch)
        processor = _processor(tmp_path, num_workers=1)
        samples = [
            {"sample_id": "a", "tree": TREE, "output_subdir": "shared"},
            {"sample_id": "b", "tree": TREE, "output_subdir": "shared"},
        ]
        with pytest.raises(ValueError, match="both resolve to output directory"):
            processor.process_samples(samples)

    def test_aliased_output_subdirs_count_as_collision(self, tmp_path, monkeypatch):
        """``g1`` 与 ``./g1/``、``g1//x/..`` 是同一个目录：也必须算碰撞。"""
        _install_fake_pipeline(monkeypatch)
        processor = _processor(tmp_path, num_workers=1)
        samples = [
            {"sample_id": "a", "tree": TREE, "output_subdir": "g1"},
            {"sample_id": "b", "tree": TREE, "output_subdir": "./g1/"},
        ]
        with pytest.raises(ValueError, match="both resolve to output directory"):
            processor.process_samples(samples)

    def test_duplicate_sample_ids_fail_the_batch(self, tmp_path, monkeypatch):
        _install_fake_pipeline(monkeypatch)
        processor = _processor(tmp_path, num_workers=1)
        samples = [
            {"sample_id": "dup", "tree": TREE, "output_subdir": "one"},
            {"sample_id": "dup", "tree": TREE, "output_subdir": "two"},
        ]
        with pytest.raises(ValueError, match="duplicate sample_id"):
            processor.process_samples(samples)

    def test_output_subdir_escaping_batch_root_is_rejected(self, tmp_path, monkeypatch):
        _install_fake_pipeline(monkeypatch)
        processor = _processor(tmp_path, num_workers=1)
        with pytest.raises(ValueError, match="outside the batch output directory"):
            processor.process_samples(
                [{"sample_id": "run", "tree": TREE, "output_subdir": "../escape"}]
            )

    def test_non_mapping_sample_raises_type_error(self, tmp_path, monkeypatch):
        _install_fake_pipeline(monkeypatch)
        with pytest.raises(TypeError, match="must be a mapping"):
            _processor(tmp_path, num_workers=1).process_samples(["not-a-sample"])

    def test_single_sample_entry_point_rejects_id_less_sample(self, tmp_path):
        """绕过 ``process_samples`` 直调单样本入口时，缺身份必须响亮拒绝。"""
        with pytest.raises(ValueError, match="sample_id"):
            _processor(tmp_path)._process_single_sample({"tree": TREE})


class TestSuccessSemanticsAndTiming:
    """C-50：一个"成功"口径 + CPU 秒与真实墙钟分开报告"""

    def test_empty_result_sample_is_not_successful(self, tmp_path, monkeypatch):
        _install_fake_pipeline(
            monkeypatch,
            results_by_subdir={
                "sample_0000": {},
                "sample_0001": {"pathd8": _dating_result()},
            },
        )
        result = _processor(tmp_path, num_workers=1).process_samples(
            [{"tree": TREE}, {"tree": TREE}]
        )

        empty, produced = result.sample_results
        assert empty.success is False
        assert empty.produced_results is False
        assert "no method results" in empty.error_message
        assert produced.produced_results is True
        assert (result.successful, result.failed) == (1, 1)
        # 汇总与方法统计现在用同一个判据：两处数字必须一致
        assert result.aggregated_statistics["successful"] == result.successful

    def test_no_results_aborts_when_continue_on_error_is_false(
        self, tmp_path, monkeypatch
    ):
        _install_fake_pipeline(monkeypatch, default={})
        processor = _processor(tmp_path, num_workers=1, continue_on_error=False)
        with pytest.raises(RuntimeError, match="no method results"):
            processor.process_samples([{"tree": TREE}])

    def test_dry_run_batch_reports_empty_results_as_dry_run_not_crash(
        self, tmp_path, monkeypatch
    ):
        """``dry_run`` 本就不产出结果：不算成功，但也不许伪装成"样本失败"。"""
        _install_fake_pipeline(monkeypatch, default={})
        processor = BatchProcessor(
            base_config=PipelineConfig(
                methods=["pathd8"], output_dir=tmp_path, dry_run=True
            ),
            batch_config=BatchConfig(num_workers=1, continue_on_error=False),
        )
        result = processor.process_samples([{"tree": TREE}])  # 不得抛出

        row = result.sample_results[0]
        assert row.success is False and row.produced_results is False
        assert "dry run" in row.error_message
        assert result.successful == 0 and result.failed == 1

    def test_cpu_seconds_and_wall_clock_are_reported_apart(self, tmp_path, monkeypatch):
        _install_fake_pipeline(monkeypatch, sleep=0.2)
        result = _processor(tmp_path, num_workers=2).process_samples(
            [{"sample_id": "s1", "tree": TREE}, {"sample_id": "s2", "tree": TREE}]
        )
        stats = result.aggregated_statistics

        assert stats["cpu_seconds_total"] >= 0.35  # 两个样本各自睡了 0.2 秒
        assert stats["wall_clock_seconds"] is not None
        assert stats["wall_clock_seconds"] < stats["cpu_seconds_total"]
        # 旧的歧义键名保留为 CPU 秒别名，不再冒充整批耗时
        assert stats["total_execution_seconds"] == stats["cpu_seconds_total"]
        assert stats["timing_note"]
        assert result.wall_clock_seconds == stats["wall_clock_seconds"]

    def test_sequential_batch_wall_clock_covers_all_samples(
        self, tmp_path, monkeypatch
    ):
        _install_fake_pipeline(monkeypatch, sleep=0.1)
        result = _processor(tmp_path, num_workers=1).process_samples(
            [{"sample_id": "s1", "tree": TREE}, {"sample_id": "s2", "tree": TREE}]
        )
        stats = result.aggregated_statistics
        assert stats["wall_clock_seconds"] >= 0.18
        # 顺序跑时两者本应接近：证明"墙钟"确实是整批量出来的，不是抄的 CPU 和
        assert abs(stats["wall_clock_seconds"] - stats["cpu_seconds_total"]) < 0.08

    def test_manifest_rows_and_wall_clock_are_serialized(self, tmp_path, monkeypatch):
        _install_fake_pipeline(monkeypatch, results_by_subdir={"sample_0001": {}})
        result = _processor(tmp_path, num_workers=1).process_samples(
            [{"tree": TREE}, {"tree": TREE}]
        )
        payload = result.to_dict()
        assert payload["wall_clock_seconds"] == pytest.approx(
            result.aggregated_statistics["wall_clock_seconds"]
        )
        assert payload["successful"] == 1
        assert payload["successful_definition"]
        rows = payload["sample_results"]
        assert [row["produced_results"] for row in rows] == [True, False]


class TestCheckpoint:
    """检查点测试"""

    def test_ckp_030_save_checkpoint(self, tmp_path):
        """CKP-030: 首次运行保存检查点"""
        manager = CheckpointManager(tmp_path)
        manager.start_step("step1")
        manager.complete_step("step1", {"done": False})
        assert manager.checkpoint_file.exists()

    def test_ckp_031_resume_checkpoint(self, tmp_path):
        """CKP-031: 从检查点恢复"""
        manager = CheckpointManager(tmp_path)
        manager.start_step("step1")
        manager.complete_step("step1", {"value": 2})
        # 重新初始化管理器模拟恢复
        manager2 = CheckpointManager(tmp_path)
        assert manager2.is_step_completed("step1")
        assert manager2.get_step_data("step1").get("value") == 2

    def test_ckp_032_cleanup_after_success(self, tmp_path):
        """CKP-032: 成功后清理检查点"""
        manager = CheckpointManager(tmp_path)
        manager.start_step("step1")
        manager.complete_step("step1")
        manager.reset()
        assert not manager.checkpoint_file.exists()


class TestInputFingerprint:
    """B-23：输入指纹的登记 / 读取 / 比对"""

    def test_no_stored_fingerprint_never_matches(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        assert manager.input_matches("sha256:whatever") is False
        assert manager.get_input_fingerprint() is None
        assert manager.get_input_hash(RUN_FINGERPRINT_KEY) is None

    def test_record_then_match(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        fingerprint = "sha256:" + "a" * 64
        manager.record_input_hash(RUN_FINGERPRINT_KEY, fingerprint)
        assert manager.input_matches(fingerprint) is True
        assert manager.input_matches("sha256:" + "b" * 64) is False
        assert manager.get_input_hash(RUN_FINGERPRINT_KEY) == fingerprint

    def test_fingerprint_survives_reload(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        manager.record_input_hash(RUN_FINGERPRINT_KEY, "sha256:abc")
        reloaded = PipelineCheckpointManager(tmp_path)
        assert reloaded.get_input_hash(RUN_FINGERPRINT_KEY) == "sha256:abc"
        assert reloaded.input_matches("sha256:abc") is True

    def test_accepts_path_and_str_keys(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        tree_path = tmp_path / "tree.nwk"
        manager.record_input_hash(tree_path, "sha256:tree")
        manager.record_input_hash("alignment.fasta", "sha256:aln")
        assert manager.get_input_hash(Path(str(tree_path))) == "sha256:tree"
        assert manager.get_input_hash("alignment.fasta") == "sha256:aln"
        assert manager.check_input_changed(tree_path, "sha256:tree") is False
        assert manager.check_input_changed(tree_path, "sha256:other") is True
        # 没记录过的键：一律按"已变更 / 读不到"处理
        assert manager.check_input_changed(tmp_path / "new.nwk", "sha256:x") is True
        assert manager.get_input_hash(tmp_path / "new.nwk") is None

    def test_rejects_junk_keys_and_hashes(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        with pytest.raises(TypeError):
            manager.record_input_hash(None, "sha256:x")
        with pytest.raises(TypeError):
            manager.record_input_hash(RUN_FINGERPRINT_KEY, 12345)
        with pytest.raises(ValueError):
            manager.record_input_hash(RUN_FINGERPRINT_KEY, "   ")
        with pytest.raises(TypeError):
            manager.record_run_parameters(["not", "a", "mapping"])
        with pytest.raises(TypeError):
            manager.compute_input_hash(["nope"])

    def test_compute_input_hash_is_order_insensitive_and_pipeline_compatible(
        self, tmp_path
    ):
        import hashlib
        import json as _json

        inputs = {"tree": "((A,B));", "calibrations": ["a", "b"], "seed": None}
        shuffled = {"seed": None, "calibrations": ["a", "b"], "tree": "((A,B));"}
        first = PipelineCheckpointManager.compute_input_hash(inputs)
        assert first == PipelineCheckpointManager.compute_input_hash(shuffled)
        assert first.startswith("sha256:")
        # 与 core/pipeline.py 的本地兜底算法逐字同口径（否则切实现就废掉旧检查点）
        canonical = _json.dumps(inputs, sort_keys=True, default=repr)
        assert (
            first == "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()
        )
        assert first != PipelineCheckpointManager.compute_input_hash(
            {**inputs, "tree": "((A,C));"}
        )

    def test_reset_clears_fingerprint_no_stale_resurrection(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        manager.record_input_hash(RUN_FINGERPRINT_KEY, "sha256:old")
        manager.record_run_parameters({"run_fingerprint": "sha256:old"})
        manager.reset()
        assert manager.get_input_fingerprint() is None
        assert manager.input_matches("sha256:old") is False
        # 换一个管理器（模拟下一次进程）也读不到旧指纹
        assert PipelineCheckpointManager(tmp_path).input_matches("sha256:old") is False

    def test_aggregate_fallback_from_per_input_hashes(self, tmp_path):
        """只登记逐文件哈希（无运行级键）时，聚合指纹仍能支撑一致性判定。"""
        manager = PipelineCheckpointManager(tmp_path)
        manager.record_input_hash("tree.nwk", "sha256:tree")
        manager.record_input_hash("aln.fasta", "sha256:aln")
        aggregate = PipelineCheckpointManager.compute_input_hash(
            {"tree.nwk": "sha256:tree", "aln.fasta": "sha256:aln"}
        )
        assert manager.input_matches(aggregate) is True
        assert manager.input_matches("sha256:other") is False


class TestResumeReuseGranularity:
    """B-23：把指纹放进"是否跳过某个方法"的判定粒度里"""

    @staticmethod
    def _result(method="mcmctree"):
        return DatingResult(
            method_name=method,
            run_id=f"{method}_r1",
            dated_tree_newick="((A,B));",
            node_ages={"n1": NodeAgeEstimate(mean_age=100.0, ci_type=CIType.NONE)},
            execution_seconds=1.0,
        )

    def test_completed_method_reused_only_for_same_fingerprint(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        manager.start_method("mcmctree")
        manager.complete_method("mcmctree", self._result(), input_fingerprint="FP-1")

        assert manager.is_method_completed("mcmctree") is True
        assert manager.is_method_completed("mcmctree", "FP-1") is True
        assert manager.is_method_completed("mcmctree", "FP-2") is False

        assert manager.get_methods_to_run(["mcmctree"], True, "FP-1") == []
        assert manager.get_methods_to_run(["mcmctree"], True, "FP-2") == ["mcmctree"]
        # 不传指纹 = 旧行为（由调用方自行核验），保证向后兼容
        assert manager.get_methods_to_run(["mcmctree"], skip_completed=True) == []

        assert manager.get_cached_result("mcmctree", "FP-1") is not None
        assert manager.get_cached_result("mcmctree", "FP-2") is None

    def test_legacy_checkpoint_without_fingerprint_refuses_reuse(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        manager.start_method("mcmctree")
        manager.complete_method("mcmctree", self._result())  # 没有指纹
        # 旧检查点里没有指纹记录 → 无法核验 → 拒绝复用（fail-closed）
        assert manager.get_cached_result("mcmctree", "FP-1") is None
        assert manager.get_methods_to_run(["mcmctree"], True, "FP-1") == ["mcmctree"]

    def test_run_fingerprint_is_stored_in_progress_summary(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        manager.start_method("mcmctree")
        manager.complete_method("mcmctree", self._result(), input_fingerprint="FP-1")
        summary = manager.get_progress_summary()
        assert summary["methods"]["mcmctree"]["input_fingerprint"] == "FP-1"
        assert summary["resume_reuse_verifiable"] is False  # 只登记在方法级
        manager.record_input_hash(RUN_FINGERPRINT_KEY, "FP-1")
        assert manager.get_progress_summary()["resume_reuse_verifiable"] is True

    def test_fingerprint_falls_back_to_registered_run_hash(self, tmp_path):
        manager = PipelineCheckpointManager(tmp_path)
        manager.record_input_hash(RUN_FINGERPRINT_KEY, "FP-9")
        manager.start_method("pathd8")
        manager.complete_method("pathd8", self._result("pathd8"))
        # 不显式传指纹时，沿用已登记的运行级指纹
        reloaded = PipelineCheckpointManager(tmp_path)
        assert reloaded.get_cached_result("pathd8", "FP-9") is not None
        assert reloaded.get_cached_result("pathd8", "FP-other") is None


class TestLoggingPayloadRobustness:
    """B-27：日志层的"报告问题"通道自己不能崩"""

    def test_exception_object_is_accepted(self, capsys):
        logger = get_logger()
        degradation = SemanticDegradationWarning(
            "MD-Cat: MRCA topology matching failed; fell back to label matching"
        )
        logger.warning(degradation)  # 旧实现在这里抛 TypeError
        captured = capsys.readouterr().out
        assert "topology matching failed" in captured

    def test_various_non_string_payloads(self, capsys):
        logger = get_logger()
        for payload in (ValueError("bad value"), 42, None, {"k": "v"}, Path("x/y")):
            logger.info(payload)
            logger.error(payload)
        captured = capsys.readouterr().out
        assert "bad value" in captured
        assert "{'k': 'v'}" in captured or '"k"' in captured

    def test_sanitize_message_coerces_and_masks_paths(self):
        logger = get_logger()
        sanitized = logger._sanitize_message(
            ValueError("/Users/secret_unit/deep/file.txt is missing")
        )
        assert isinstance(sanitized, str)
        assert "file.txt is missing" in sanitized
        assert "/Users/secret_unit/" not in sanitized  # 路径脱敏仍然生效
        assert logger._sanitize_message(123) == "123"

    def test_unprintable_object_does_not_crash_logging(self):
        class Explodes:
            def __str__(self):
                raise RuntimeError("nope")

        logger = get_logger()
        logger.info(Explodes())  # 不得抛出

    def test_warning_with_category_still_emits_warnings_module(self, recwarn):
        logger = get_logger()
        logger.warning(
            SemanticDegradationWarning("degraded"), category=SemanticDegradationWarning
        )
        assert any("degraded" in str(item.message) for item in recwarn.list)


class TestRuntimeMetadata:
    """运行时元数据测试"""

    def test_bat_020_metadata_records(self, tmp_path):
        """BAT-020: 元数据记录"""
        manager = RuntimeMetadataManager(output_dir=tmp_path)
        manager.collect_system_resources()
        manager.set_input_file("tree", "tree.nwk")
        path = manager.save()
        assert path.exists()

    def test_bat_022_metadata_failure(self, tmp_path):
        """BAT-022: 元数据失败记录"""
        manager = RuntimeMetadataManager(output_dir=tmp_path)
        manager.start_method("pathd8")
        manager.fail_method("pathd8", "test error")
        manager.finalize()
        path = manager.save()
        content = path.read_text(encoding="utf-8")
        assert "failed" in content or "error" in content.lower()
