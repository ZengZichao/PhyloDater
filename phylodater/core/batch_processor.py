"""
Batch Processing - 批量处理接口

支持多样本、多基因的批量定年分析。
设计用于处理大规模物种水平的时间树构建。
"""

import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from ..infrastructure import get_logger
from ..infrastructure.safe_io import safe_writer
from ..models import DatingResult, PhylogeneticTree
from .pipeline import DatingPipeline, PipelineConfig


@dataclass
class BatchConfig:
    """批量处理配置"""

    num_workers: int = 4
    continue_on_error: bool = True
    collect_summary_statistics: bool = True
    output_aggregated_results: bool = True


@dataclass
class SampleResult:
    """单个样本的定年结果

    ``output_dir`` 指向该样本自己的输出目录——真正的年龄数据只在其中的
    ``comparison_table.tsv`` 里，本清单（``BatchResult.to_dict``）不含年龄。
    """

    sample_id: str
    success: bool
    results: Optional[Dict[str, DatingResult]] = None
    error_message: Optional[str] = None
    execution_seconds: float = 0.0
    output_dir: Optional[Path] = None

    @property
    def produced_results(self) -> bool:
        """C-50：批量"成功"的**唯一**口径——跑通且真的产出了方法结果。

        旧实现里 ``process_samples`` 用 ``if r.success`` 计数、``_compute_statistics``
        用 ``if r.success and r.results`` 过滤，同一个词在两处指两件事：汇总印
        ``Successful: 40``、方法统计却只有 ``34 samples``，差的那 6 个样本其实什么
        都没产出（``DatingPipeline.run`` 在方法全部失败时未必抛出）。判据收拢到这一
        个属性上，两处共用，口径不可能再分叉。
        """
        return bool(self.success and self.results)


#: ``batch_results.json`` 顶层的显式声明（审阅项 C-51）：这份清单**不含年龄**。
MANIFEST_NOTE = (
    "This manifest contains per-sample status and node-count statistics ONLY. "
    "It has no age estimates: node ages live in each sample's own "
    "<output_dir>/comparison_table.tsv (topology-aligned, with per-method "
    "*_ci_type columns) and <output_dir>/*.nhx dated trees."
)


@dataclass
class BatchResult:
    """批量处理结果

    注意：``to_dict()`` 故意不序列化 ``results``（真正的 ``DatingResult`` 对象），
    因此写出的 ``batch_results.json`` 是一份**清单/统计**，不是年龄汇总（C-51）。
    顶层的 ``note`` / ``age_data_location`` 与每个样本的 ``output_dir`` 一起，
    把下游读者指向真正有年龄的文件。
    """

    total_samples: int
    successful: int
    failed: int
    sample_results: List[SampleResult]
    aggregated_statistics: Optional[Dict[str, Any]] = None
    #: C-50：整批的**真实墙钟**（秒）。各样本墙钟之和是 CPU 秒语义，两者在
    #: ``num_workers > 1`` 时相差可达 worker 数量级，必须分开报告。
    wall_clock_seconds: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "artifact_type": "batch_manifest",
            "note": MANIFEST_NOTE,
            "age_data_location": "<sample_results[].output_dir>/comparison_table.tsv",
            "total_samples": self.total_samples,
            "successful": self.successful,
            "failed": self.failed,
            "successful_definition": "sample ran AND produced at least one method result",
            "wall_clock_seconds": self.wall_clock_seconds,
            "sample_results": [
                {
                    "sample_id": sr.sample_id,
                    "success": sr.success,
                    "produced_results": sr.produced_results,
                    "error_message": sr.error_message,
                    "execution_seconds": sr.execution_seconds,
                    "output_dir": str(sr.output_dir) if sr.output_dir else None,
                    "age_data_file": (
                        str(sr.output_dir / "comparison_table.tsv")
                        if sr.output_dir
                        else None
                    ),
                    "methods": sorted(sr.results.keys()) if sr.results else [],
                    "node_key_counts": (
                        {
                            method: len(getattr(result, "node_ages", {}) or {})
                            for method, result in sr.results.items()
                        }
                        if sr.results
                        else {}
                    ),
                }
                for sr in self.sample_results
            ],
            "aggregated_statistics": self.aggregated_statistics,
        }


class BatchProcessor:
    """
    批量处理器

    支持多样本并行定年分析，适用于：
    - 多个基因的独立定年
    - 多个物种树的批量分析
    - 超大规模时间树构建

    设计原则：
    - 每个样本独立运行完整管线
    - 每个样本有**互不相同**的身份：``sample_id`` 与输出目录（含其下的检查点文件）
      一起构成续跑判定的键，任一碰撞都可能把 A 样本的缓存结果交给 B 样本（C-49）
    - 结果缓存避免重复计算
    - 支持断点续跑
    """

    def __init__(
        self, base_config: PipelineConfig, batch_config: Optional[BatchConfig] = None
    ) -> None:
        self.base_config = base_config
        self.batch_config = batch_config or BatchConfig()
        self.logger = get_logger()

    def process_samples(
        self,
        samples: List[Dict[str, Any]],
        progress_callback: Optional[Callable[[int, int, str], None]] = None,
    ) -> BatchResult:
        """
        批量处理样本

        Args:
            samples: 样本列表，每个样本包含：
                - sample_id: 样本唯一标识（缺失时按索引派生 ``sample_{i:04d}``，见 C-49）
                - tree: PhylogeneticTree 对象或 Newick 字符串
                - alignment_path: 比对文件路径
                - calibrations: 校准点列表
                - output_subdir: 可选的输出子目录（整批必须互不相同）
            progress_callback: 进度回调函数 (completed, total, current_sample_id)

        Returns:
            BatchResult: 批量处理结果

        Raises:
            ValueError: 两个样本落到同一个输出目录（或共用同一个 ``sample_id``）
        """
        # C-49：身份必须在分发**之前**敲定。缺 sample_id 的旧写法
        # （``sample.get("sample_id", "unknown")``）会让多个样本共用同一输出目录
        # **和同一个 checkpoint 文件**（``_create_sample_config`` 无条件开启检查点，
        # 而检查点就建在 output_dir 里）；叠上 B-23 的"按输入指纹续跑"，第二个样本
        # 会把第一个样本的缓存结果当成自己的结果并双双报成功。
        samples = self._resolve_sample_identities(samples)
        total = len(samples)
        sample_results: List[SampleResult] = []

        self.logger.section(f"Batch Processing: {total} samples")

        # C-50：真实墙钟只在这里量——各样本墙钟之和是 CPU 秒，并行时可达 num_workers 倍。
        wall_start = time.perf_counter()
        if self.batch_config.num_workers > 1:
            sample_results = self._process_parallel(samples, progress_callback)
        else:
            sample_results = self._process_sequential(samples, progress_callback)
        wall_clock_seconds = time.perf_counter() - wall_start

        successful = sum(1 for r in sample_results if r.produced_results)
        failed = total - successful

        batch_result = BatchResult(
            total_samples=total,
            successful=successful,
            failed=failed,
            sample_results=sample_results,
            wall_clock_seconds=wall_clock_seconds,
        )

        if self.batch_config.collect_summary_statistics:
            batch_result.aggregated_statistics = self._compute_statistics(
                sample_results, wall_clock_seconds=wall_clock_seconds
            )

        self._log_summary(batch_result)

        return batch_result

    # ------------------------------------------------------------------ #
    # C-49：样本身份与输出目录唯一性
    # ------------------------------------------------------------------ #

    def _resolve_sample_identities(
        self, samples: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """给每个样本敲定唯一 ``sample_id`` 与唯一 ``output_subdir``，否则立即失败。

        返回的是**浅拷贝**列表（不修改调用方传入的字典），其中 ``sample_id`` 与
        ``output_subdir`` 都已归一化，后续执行链因此只有一个身份来源。
        """
        resolved: List[Dict[str, Any]] = []
        seen_ids: Dict[str, int] = {}
        seen_dirs: Dict[Path, Tuple[int, str]] = {}

        for index, sample in enumerate(samples):
            if not isinstance(sample, dict):
                raise TypeError(
                    f"Batch sample #{index} must be a mapping (with a 'sample_id' key), "
                    f"got {type(sample).__name__}"
                )

            entry = dict(sample)
            sample_id = self._normalize_sample_id(entry.get("sample_id"), index)
            if sample_id in seen_ids:
                raise ValueError(
                    f"C-49: duplicate sample_id {sample_id!r} at batch indices "
                    f"{seen_ids[sample_id]} and {index}. Sample ids are the checkpoint "
                    "identity, so two samples sharing one id can be resumed into each "
                    "other's cached results."
                )
            seen_ids[sample_id] = index
            entry["sample_id"] = sample_id

            output_subdir = self._normalize_output_subdir(
                entry.get("output_subdir"), sample_id
            )
            output_dir = self._resolve_output_dir(output_subdir, sample_id)
            if output_dir in seen_dirs:
                prev_index, prev_id = seen_dirs[output_dir]
                raise ValueError(
                    f"C-49: batch samples {prev_id!r} (index {prev_index}) and "
                    f"{sample_id!r} (index {index}) both resolve to output directory "
                    f"{output_dir}. They would overwrite each other's artifacts and "
                    "share one checkpoint file — give them distinct 'output_subdir' "
                    "values (or distinct sample_ids)."
                )
            seen_dirs[output_dir] = (index, sample_id)
            entry["output_subdir"] = output_subdir

            resolved.append(entry)

        return resolved

    def _normalize_sample_id(self, raw: Any, index: int) -> str:
        """样本身份：有则用、空则按索引派生（并在清单里留下可追溯的名字）。"""
        sample_id = "" if raw is None else str(raw).strip()
        if not sample_id:
            sample_id = f"sample_{index:04d}"
            self.logger.warning(
                f"C-49: batch sample #{index} has no usable 'sample_id'; derived "
                f"{sample_id!r} from its index so it cannot share another sample's "
                "output directory or checkpoint."
            )
        return sample_id

    @staticmethod
    def _normalize_output_subdir(raw: Any, sample_id: str) -> str:
        subdir = "" if raw is None else str(raw).strip()
        return subdir or sample_id

    def _resolve_output_dir(self, output_subdir: str, sample_id: str) -> Path:
        """把 ``output_subdir`` 归一化成基准目录下的绝对路径（唯一性判据）。

        用解析后的路径而不是原始字符串比较，才能让 ``a/b`` 与 ``a//b``、``./a`` 这类
        同义写法也算碰撞；顺带拦住跑到基准目录之外（``..``／绝对路径）的样本，否则
        "唯一"只是字符串层面的假象。
        """
        base = Path(self.base_config.output_dir).expanduser().resolve()
        candidate = (base / output_subdir).resolve()
        if candidate != base and base not in candidate.parents:
            raise ValueError(
                f"C-49: sample {sample_id!r} resolves outside the batch output "
                f"directory ({candidate} not under {base}); batch outputs must stay "
                "isolated per sample."
            )
        return candidate

    @staticmethod
    def _require_sample_id(sample: Dict[str, Any]) -> str:
        """单样本入口的身份守卫：``process_samples`` 已派生好，直接调用必须自带。"""
        sample_id = str(sample.get("sample_id") or "").strip()
        if not sample_id:
            raise ValueError(
                "C-49: batch samples must carry a non-empty 'sample_id' — it is both "
                "the output directory and the checkpoint identity. Run the batch "
                "through process_samples() (which derives sample_{i:04d}) instead of "
                "calling _process_single_sample() with an id-less sample."
            )
        return sample_id

    def _process_sequential(
        self,
        samples: List[Dict[str, Any]],
        progress_callback: Optional[Callable[[int, int, str], None]],
    ) -> List[SampleResult]:
        """顺序处理样本"""
        results = []
        for i, sample in enumerate(samples):
            result = self._process_single_sample(sample)
            results.append(result)

            if progress_callback:
                progress_callback(i + 1, len(samples), sample["sample_id"])

        return results

    def _process_parallel(
        self,
        samples: List[Dict[str, Any]],
        progress_callback: Optional[Callable[[int, int, str], None]],
    ) -> List[SampleResult]:
        """并行处理样本"""
        results: Dict[int, SampleResult] = {}
        lock = threading.Lock()
        completed_count = [0]

        def process_with_index(
            args: Tuple[int, Dict[str, Any]],
        ) -> SampleResult:
            i, sample = args
            result = self._process_single_sample(sample)

            with lock:
                results[i] = result
                completed_count[0] += 1
                if progress_callback:
                    progress_callback(
                        completed_count[0], len(samples), sample["sample_id"]
                    )

            return result

        with ThreadPoolExecutor(max_workers=self.batch_config.num_workers) as executor:
            futures = {
                executor.submit(process_with_index, (i, sample)): i
                for i, sample in enumerate(samples)
            }

            for future in as_completed(futures):
                future.result()

        return [results[i] for i in sorted(results.keys())]

    def _process_single_sample(self, sample: Dict[str, Any]) -> SampleResult:
        """处理单个样本"""
        sample_id = self._require_sample_id(sample)
        start_time = time.perf_counter()
        sample_output_dir: Optional[Path] = None

        try:
            output_subdir = self._normalize_output_subdir(
                sample.get("output_subdir"), sample_id
            )
            config = self._create_sample_config(output_subdir)
            sample_output_dir = config.output_dir

            tree = sample.get("tree")
            if isinstance(tree, str):
                tree = PhylogeneticTree.from_newick(tree)
            elif isinstance(tree, Path):
                tree = PhylogeneticTree.from_file(tree)
            if not isinstance(tree, PhylogeneticTree):
                # 批量入口的样本来自 YAML/JSON，"没写 tree" 或"写了个不能识别的
                # 形状"都很常见：当场报清楚是哪个样本，而不是把 None 递进管线，
                # 由深处的属性访问报一个看不出来源的 AttributeError。
                raise TypeError(
                    f"sample '{sample_id}': 'tree' must be a PhylogeneticTree, a "
                    f"Newick string or a path; got {type(tree).__name__}: {tree!r}"
                )

            alignment_path = sample.get("alignment_path")
            calibrations = sample.get("calibrations", [])

            pipeline = DatingPipeline(config)
            results = pipeline.run(tree, alignment_path, calibrations)

            if not results:
                # C-50：管线"跑完了"却不产出任何方法结果时未必抛出（C-8 的那批
                # ``except Exception``），旧口径把这种样本计成成功。这里把它当成失败
                # 抛出，于是 ``success=True`` 与"有结果"在结构上等价，且
                # ``continue_on_error=False`` 也能真的中止。
                if config.dry_run:
                    # dry-run 依设计不产出结果：既不该计成"成功产出了东西"，也不该
                    # 被记成样本失败——单列一条自解释的结果（清单里看得见原因）。
                    return SampleResult(
                        sample_id=sample_id,
                        success=False,
                        error_message=(
                            "dry run: the pipeline produces no results by design "
                            "(PipelineConfig.dry_run=True)"
                        ),
                        execution_seconds=time.perf_counter() - start_time,
                        output_dir=sample_output_dir,
                    )
                raise RuntimeError(
                    "Dating pipeline produced no method results "
                    "(every requested method failed or was skipped)"
                )

            execution_seconds = time.perf_counter() - start_time

            return SampleResult(
                sample_id=sample_id,
                success=True,
                results=results,
                execution_seconds=execution_seconds,
                output_dir=sample_output_dir,
            )

        except Exception as e:
            execution_seconds = time.perf_counter() - start_time
            self.logger.error(f"Sample {sample_id} failed: {e}")

            if self.batch_config.continue_on_error:
                return SampleResult(
                    sample_id=sample_id,
                    success=False,
                    error_message=str(e),
                    execution_seconds=execution_seconds,
                    output_dir=sample_output_dir,
                )
            else:
                raise

    def _create_sample_config(self, output_subdir: str) -> PipelineConfig:
        """为单个样本创建配置"""
        import copy

        config = copy.deepcopy(self.base_config)
        config.output_dir = self.base_config.output_dir / output_subdir
        config.enable_checkpoint = True

        return config

    def _compute_statistics(
        self,
        sample_results: List[SampleResult],
        wall_clock_seconds: Optional[float] = None,
    ) -> Dict[str, Any]:
        """计算汇总统计

        ``wall_clock_seconds`` 是整批的真实耗时（由 :meth:`process_samples` 在分发
        外层量得）。各样本墙钟之和是 **CPU 秒**语义：``num_workers=4`` 时它约为真实
        耗时的 4 倍，所以两个数必须分开报告（C-50）。
        """
        successful_results = [r for r in sample_results if r.produced_results]

        if not successful_results:
            return {"error": "No successful results to aggregate"}

        all_methods: Set[str] = set()
        for result in successful_results:
            if result.results:
                all_methods.update(result.results.keys())

        method_stats = {}
        for method in all_methods:
            method_results = [
                r.results[method]
                for r in successful_results
                if r.results and method in r.results
            ]

            if method_results:
                node_counts = [len(r.node_ages) for r in method_results]
                method_stats[method] = {
                    "count": len(method_results),
                    "avg_nodes": (
                        sum(node_counts) / len(node_counts) if node_counts else 0
                    ),
                    "min_nodes": min(node_counts) if node_counts else 0,
                    "max_nodes": max(node_counts) if node_counts else 0,
                }

        cpu_seconds_total = sum(r.execution_seconds for r in sample_results)
        stats: Dict[str, Any] = {
            "total_samples": len(sample_results),
            # C-50：与 ``BatchResult.successful`` 同一个判据（产出非空结果）。
            "successful": len(successful_results),
            "methods": method_stats,
            "cpu_seconds_total": cpu_seconds_total,
            # 旧键名保留为 CPU 秒别名，避免打断已有读者，同时由 ``timing_note`` 说明语义。
            "total_execution_seconds": cpu_seconds_total,
            "avg_cpu_seconds_per_sample": (
                cpu_seconds_total / len(sample_results) if sample_results else 0.0
            ),
            "avg_execution_seconds": (
                cpu_seconds_total / len(sample_results) if sample_results else 0.0
            ),
            "wall_clock_seconds": wall_clock_seconds,
            "num_workers": self.batch_config.num_workers,
            "timing_note": (
                "cpu_seconds_total is the sum of per-sample wall-clocks and therefore "
                "CPU-second semantics (up to ~num_workers x larger than elapsed time "
                "under parallel workers); wall_clock_seconds is the real elapsed time "
                "of the whole batch."
            ),
        }
        if wall_clock_seconds and wall_clock_seconds > 0:
            stats["parallel_efficiency"] = cpu_seconds_total / (
                wall_clock_seconds * max(1, self.batch_config.num_workers)
            )
        return stats

    def _log_summary(self, result: BatchResult) -> None:
        """记录汇总信息"""
        self.logger.section("Batch Processing Summary")
        self.logger.kv("Total Samples", result.total_samples)
        self.logger.kv("Successful (produced results)", result.successful)
        self.logger.kv("Failed", result.failed)

        if result.aggregated_statistics:
            stats = result.aggregated_statistics
            if "methods" in stats:
                self.logger.info("Method Statistics:")
                for method, method_stats in stats["methods"].items():
                    self.logger.kv(
                        f"  {method}",
                        f"{method_stats['count']} samples, "
                        f"avg nodes: {method_stats['avg_nodes']:.1f}",
                    )
            if "cpu_seconds_total" in stats:
                self.logger.kv(
                    "CPU seconds (sum of per-sample)",
                    f"{stats['cpu_seconds_total']:.1f}",
                )
            if stats.get("wall_clock_seconds") is not None:
                self.logger.kv(
                    "Wall clock (whole batch)", f"{stats['wall_clock_seconds']:.1f}"
                )
            elif result.wall_clock_seconds is not None:
                self.logger.kv(
                    "Wall clock (whole batch)", f"{result.wall_clock_seconds:.1f}"
                )

    def save_results(self, result: BatchResult, output_path: Path) -> None:
        """保存批量处理**清单**（状态 + 统计，不含年龄数据）。

        C-51：``results``（真正的 ``DatingResult``）有意不序列化，所以这个 JSON
        里没有任何节点年龄。清单顶层写了 ``artifact_type`` / ``note``，逐样本写了
        ``output_dir`` 与 ``age_data_file``，避免下游把 ``batch_results.json`` 当成
        年龄汇总表去算平均值。若希望产物名本身就说明这一点，请把 ``output_path``
        取为 ``batch_manifest.json``。
        """
        import json

        output_path.parent.mkdir(parents=True, exist_ok=True)

        with safe_writer(output_path, encoding="utf-8") as f:
            json.dump(result.to_dict(), f, indent=2, ensure_ascii=False)

        self.logger.info(f"Batch results saved to {output_path}")
        self.logger.info(
            "Note: this manifest holds status/statistics only — node ages are in each "
            "sample's <output_dir>/comparison_table.tsv (see the 'age_data_file' field)."
        )
        if str(output_path).endswith("batch_results.json"):
            self.logger.debug(
                "C-51: 'batch_results.json' is a manifest, not an age summary; "
                "consider naming it 'batch_manifest.json'."
            )
