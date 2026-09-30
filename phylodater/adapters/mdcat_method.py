"""
MDCatMethod - MD-Cat 定年适配器

实现 MD-Cat (EM 算法分类速率模型) 定年功能
根据诊断报告：
- 使用 treeswift 库
- 默认注入 -b 参数（化石定年模式）
- 内部节点预打标签
- 支持置信区间计算 (--CI)
- 支持直接 Python 导入
"""

import re
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from ..core import DatingMethod, DatingMethodRegistry
from ..core.exceptions import (
    CalibrationError,
    ExecutionError,
    ResultParsingError,
    SemanticDegradationWarning,
)
from ..core.method_interface import InputFileInfo
from ..infrastructure import ProcessRunner, get_logger
from ..infrastructure.configuration import (
    CommonConfig,
    MDCatConfig,
    SoftwarePaths,
    ToolConfig,
)
from ..infrastructure.safe_io import safe_writer
from ..models import (
    CalibrationPoint,
    CIType,
    DatingResult,
    NodeAgeEstimate,
    PhylogeneticTree,
    strip_leading_newick_comments,
)


class MDCatMethod(DatingMethod[MDCatConfig]):
    """
    MD-Cat 定年适配器

    特性：
    - 分类速率模型 (Categorical Rate Model) / EM 算法
    - 默认 -b 参数（化石定年，时间向后）
    - 内部节点预打标签
    - 支持置信区间 (--CI)
    - 支持直接 Python 导入（更高效）
    """

    def __init__(
        self,
        config: Union[ToolConfig, MDCatConfig],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        # 入口保留双形态（历史上两种调用都存在），但在交给基类之前先归一：
        # 基类把 ``self.config`` 记为 ``MDCatConfig``，适配器后续只读子配置。
        # ``isinstance`` 与原先的 ``hasattr(config, "mdcat")`` 对这两种入参等价：
        # ``ToolConfig`` 必有 ``mdcat``，``MDCatConfig`` 必无。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            method_config = config.mdcat
            if common is None:
                common = config.common
        else:
            method_config = config
        super().__init__(
            method_config, output_dir, software_paths, common_config=common
        )

        self.logger = get_logger()
        # ``_input_files`` 同时存文件路径与执行统计量（seq_len / execution_time /
        # log_likelihood），详见 ``InputFileInfo`` 的说明。
        self._input_files: Dict[str, InputFileInfo] = {}
        self._calibrations: List[CalibrationPoint] = []
        self._node_label_map: Dict[str, str] = {}
        self._use_direct_import = self.config.use_direct_import
        self._mdcat_func: Optional[Callable[..., Any]] = None

    @property
    def method_name(self) -> str:
        return "mdcat"

    def _get_executable(self) -> str:
        """获取 MD-Cat 可执行文件路径，优先使用 software_paths 自定义路径。"""
        path = self.get_software_path("mdcat_bin")
        return path if path else self.config.mdcat_bin

    def validate_environment(self) -> bool:
        """检测 MD-Cat 是否可用"""
        # 检查直接导入可用性
        if self._use_direct_import:
            try:
                from emd.emd_normal_lib import MDCat as MDCat_func

                self._mdcat_func = MDCat_func
                self.software_version = self._package_version("MD-Cat")
                self.logger.info(
                    "MD-Cat (emd) Python package detected - will use direct import"
                )
                return True
            except ImportError:
                self.logger.warning(
                    "MD-Cat direct import not available, falling back to subprocess"
                )
                self._use_direct_import = False

        # 检查命令行可用性
        runner = ProcessRunner()
        executable = self._get_executable()
        available = runner.check_executable(executable)

        if available:
            self.software_version = self._package_version("MD-Cat") or executable
            self.logger.info("MD-Cat CLI detected")
        else:
            self.logger.warning(f"MD-Cat not found: {executable}")

        return available

    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
        生成 MD-Cat 输入文件

        1. 内部节点预打标签
        2. 生成约束文件
        3. 获取序列长度
        """
        self._calibrations = calibrations

        # 获取序列长度
        seq_len = 1000
        if alignment_path and alignment_path.exists():
            from ..infrastructure import AlignmentMetadataExtractor

            extractor = AlignmentMetadataExtractor()
            metadata = extractor.extract(alignment_path)
            seq_len = metadata.sequence_length

        # 内部节点预打标签（保持与全平台其它方法一致的 internal_node_N 命名）
        tree_with_labels = self._label_internal_nodes(tree, calibrations)
        tree_file = self.work_dir / "mdcat_tree.nwk"
        tree_with_labels.write(tree_file)

        # 生成约束文件
        constraints_file = self.work_dir / "mdcat_constraints.txt"

        written = 0
        with safe_writer(constraints_file) as f:
            for cal in calibrations:
                if not cal.resolved_taxa:
                    # C-22：被跳过的校准必须留痕（与 wlogdate 对齐）。
                    self.logger.warning(
                        f"MD-Cat: skipping calibration '{cal.name}': no resolved_taxa"
                    )
                    continue

                # 使用所有有效叶节点
                taxa_str = "+".join(cal.resolved_taxa)

                constraint = cal.age_constraint
                if constraint is None:
                    # 没有年龄界的校准点无法渲染成 MD-Cat 的 mrca 行：与上面
                    # 缺 resolved_taxa 的处理保持一致——留痕并跳过，而不是
                    # 在 ``None.to_software_format()`` 上碰一个 AttributeError。
                    # 若一条也没写出去，下面的 written==0 会硬失败。
                    self.logger.warning(
                        f"MD-Cat: skipping calibration '{cal.name}': no age_constraint"
                    )
                    continue

                # 统一使用约束对象的 to_software_format 生成 MD-Cat 格式
                formatted = constraint.to_software_format(
                    "mdcat", name=cal.name, taxa_str=taxa_str
                )
                f.write(f"{formatted}\n")
                written += 1

        # C-22：一条校准都没写出去时必须显式失败，否则 md-cat 会带着空约束文件继续跑。
        if written == 0:
            raise CalibrationError(
                "MD-Cat: no valid calibrations were written to the constraint file. "
                "Check that each calibration resolves to at least one leaf "
                "(mrca_leaf_pair / resolved_taxa)."
            )

        self._input_files = {
            "tree": tree_file,
            "constraints": constraints_file,
            "seq_len": seq_len,
        }

        self.logger.info("Generated MD-Cat input files")

        return self._input_files

    def _label_internal_nodes(
        self, tree: PhylogeneticTree, calibrations: List[CalibrationPoint]
    ) -> PhylogeneticTree:
        """为内部节点预打标签，并建立与 PhyloDater 统一命名对齐的映射。

        MD-Cat 需要所有内部节点都有唯一标签。本方法：
        - 校准节点直接以校准点名作为标签；
        - 非校准内部节点以 ``IntNode*`` 作为处理标签，但映射到
          ``internal_node_{unnamed_idx}``（从 1 开始，跳过校准节点），
          与 pathd8/r8s/treepl/mcmctree/lsd2 等方法的未命名节点命名一致。
        """
        try:
            import treeswift

            ts_tree = treeswift.read_tree_newick(tree.newick)

            # 建立节点 → 校准名的映射（根节点或 mrca_leaf_pair 的 MRCA）
            cal_name_by_node: Dict[Any, str] = {}
            for cal in calibrations:
                if not cal.name:
                    continue
                try:
                    if cal.is_root_node:
                        cal_name_by_node[ts_tree.root] = cal.name
                    elif cal.mrca_leaf_pair and len(cal.mrca_leaf_pair) == 2:
                        tip1, tip2 = cal.mrca_leaf_pair
                        leaf1 = next(
                            (n for n in ts_tree.traverse_leaves() if n.label == tip1),
                            None,
                        )
                        leaf2 = next(
                            (n for n in ts_tree.traverse_leaves() if n.label == tip2),
                            None,
                        )
                        if leaf1 is not None and leaf2 is not None:
                            mrca = ts_tree.mrca({leaf1.label, leaf2.label})
                            if mrca is not None:
                                cal_name_by_node[mrca] = cal.name
                except Exception as e:
                    # C-59：单个校准点定位失败要升级为 warning，否则它没标签、
                    # 下游按位置匹配精度下降却无人知晓。
                    self.logger.warning(
                        f"MD-Cat: could not resolve calibration node for '{cal.name}': "
                        f"{e}. This calibration has no internal-node label and may fall "
                        f"back to lower-precision position matching in the output."
                    )

            # 先按前序遍历给所有内部节点分配处理标签，同时记录映射
            processing_idx = 0
            unnamed_idx = 1
            self._node_label_map.clear()
            for node in ts_tree.traverse_preorder():
                if node.is_leaf():
                    continue
                cal_name = cal_name_by_node.get(node)
                if cal_name:
                    node.label = cal_name
                    # 校准名保持原样；chronogram 解码时也会通过 label_to_calname 命中
                    self._node_label_map[cal_name] = cal_name
                else:
                    temp_label = f"IntNode{processing_idx}"
                    node.label = temp_label
                    self._node_label_map[temp_label] = f"internal_node_{unnamed_idx}"
                    processing_idx += 1
                    unnamed_idx += 1

            new_newick = ts_tree.newick()
            return PhylogeneticTree(new_newick)

        except ImportError:
            # C-59：MD-Cat 需要内部节点唯一标签；缺 treeswift 时按位置匹配精度下降，
            # 用显式的降级警告提示，而不是静默返回未打标签的树。
            self.logger.warning(
                "treeswift not available: MD-Cat internal-node labeling skipped. The "
                "output tree will lack calibration labels and parse_results falls back "
                "to lower-precision position/chronogram matching.",
                category=SemanticDegradationWarning,
            )
            return tree

    def execute(self) -> bool:
        """执行 MD-Cat 分析"""
        if self._use_direct_import and self._mdcat_func:
            return self._execute_direct()
        else:
            return self._execute_subprocess()

    def _execute_direct(self) -> bool:
        """使用直接 Python 导入执行（更高效）"""
        self.logger.info("Running MD-Cat via direct import (more efficient)")

        # 进到这里就说明调用方要求走导入通路；函数没加载好是真实缺陷，
        # 直接抛成 "'NoneType' object is not callable" 会让用户无从下手。
        if self._mdcat_func is None:
            raise ExecutionError(
                "MD-Cat: direct-import execution requested but the MDCat callable "
                "is not loaded; call validate_environment() first, or set "
                "software.mdcat.use_direct_import=false to use the md-cat CLI."
            )

        try:
            import treeswift

            tree_file = self._input_files["tree"]
            constraints_file = self._input_files["constraints"]
            seq_len = self._require_input_int("seq_len")

            # 自适应计算量：MD-Cat EM 大致为 O(nrep * max_iter * seq_len * ncat)。
            # 默认参数（nrep=100, max_iter=100, ncat=50）在较大蛋白比对上会非常慢。
            # 这里按序列长度对 nrep / max_iter 做分段缩放；对超过 20000 位点的大
            # 数据进一步降低，避免不可接受的运行时间。
            nrep = self.config.nrep
            max_iter = self.config.max_iter
            if seq_len > 0:
                if seq_len >= 50000:
                    nrep = min(nrep, 3)
                    max_iter = min(max_iter, 5)
                elif seq_len >= 20000:
                    nrep = min(nrep, 5)
                    max_iter = min(max_iter, 10)
                elif seq_len >= 5000:
                    nrep = min(nrep, 10)
                    max_iter = min(max_iter, 20)
                # seq_len < 5000 保持配置默认值不变
                if nrep != self.config.nrep or max_iter != self.config.max_iter:
                    self.logger.info(
                        f"MD-Cat adaptive scaling for seq_len={seq_len}: "
                        f"nrep={nrep} (cfg={self.config.nrep}), "
                        f"max_iter={max_iter} (cfg={self.config.max_iter})"
                    )

            # 读取树
            ts_tree = treeswift.read_tree_newick(str(tree_file))

            # 确定时间参数
            if self.config.backward_time:
                tR = 1 if not self.config.as_date else None
                tL = 0 if not self.config.as_date else None
            else:
                tR = 0 if not self.config.as_date else None
                tL = 1 if not self.config.as_date else None

            # 解析随机种子（使用 common 配置）
            randseed = self.common_config.seed

            # CI 选项
            ci_nboots = self.config.ci_nboots
            CI_options = None
            if ci_nboots > 0:
                # bootstrap 在大数据上同样耗时，一并缩放；
                # 经验上 ≥20000 位点时即使少量 bootstrap 也可能极慢，直接关闭。
                if seq_len >= 20000:
                    ci_nboots = 0
                    self.logger.info(
                        f"MD-Cat CI disabled for seq_len={seq_len} (bootstrap too slow); "
                        "set ci_nboots manually to override"
                    )
                elif seq_len >= 5000:
                    ci_nboots = min(ci_nboots, 5)
                if ci_nboots > 0:
                    CI_options = {
                        "nboots": ci_nboots,
                        "p_lower": self.config.ci_plower,
                        "p_upper": self.config.ci_pupper,
                    }

            # 注释选项
            place_mu = self.config.annotate_level >= 2
            place_q = self.config.annotate_level >= 3

            start_time = time.time()

            # 调用 MDCat 函数
            best_tree, best_llh, best_phi, best_omega = self._mdcat_func(
                ts_tree,
                self.config.ncat,
                sampling_time=str(constraints_file),
                s=seq_len,
                nrep=nrep,
                maxIter=max_iter,
                refTree=None,
                fixed_tau=False,
                fixed_omega=False,
                verbose=self.common_config.verbose,
                pseudo=self.config.pseudocount,
                randseed=randseed,
                place_mu=place_mu,
                place_q=place_q,
                init_Q=None,
                root_time=tR,
                leaf_time=tL,
                bw_time=self.config.backward_time,
                as_date=self.config.as_date,
                CI_options=CI_options,
            )

            elapsed = time.time() - start_time

            # 写出结果树
            output_tree = self.work_dir / "mdcat_output.tre"
            best_tree.write_tree_newick(str(output_tree))

            self._input_files["output_tree"] = output_tree
            self._input_files["execution_time"] = elapsed
            self._input_files["log_likelihood"] = best_llh

            self.logger.success(
                f"MD-Cat execution completed in {elapsed:.2f}s, log-likelihood: {best_llh}"
            )
            return True

        except (ImportError, AttributeError) as e:
            # C-58：只有"库不可用/上游 API 缺属性"才回退到 CLI，并点明真实成因；
            # 其余（TypeError 签名漂移、计算错误等）上抛，避免把真实错误伪装成 CLI 失败。
            self.logger.warning(
                f"MD-Cat direct-import path unavailable ({type(e).__name__}: {e}); "
                "falling back to the md-cat CLI (parameters may be auto-scaled)."
            )
            return self._execute_subprocess()

    def _execute_subprocess(self) -> bool:
        """使用 MD-Cat CLI (md-cat) 执行

        注意：MD-Cat 的可执行文件名为 ``md-cat``，参数与 emd 包内部的 Python API
        不同。这里使用正确的 CLI 参数调用，并对大序列比对做自适应缩放，避免
        默认参数下运行时间过长。
        """
        tree_file = self._input_files["tree"]
        constraints_file = self._input_files["constraints"]
        seq_len = self._require_input_int("seq_len")

        # 自适应缩放（与 direct import 保持一致）
        nrep = self.config.nrep
        max_iter = self.config.max_iter
        if seq_len > 0:
            if seq_len >= 50000:
                nrep = min(nrep, 3)
                max_iter = min(max_iter, 5)
            elif seq_len >= 20000:
                nrep = min(nrep, 5)
                max_iter = min(max_iter, 10)
            elif seq_len >= 5000:
                nrep = min(nrep, 10)
                max_iter = min(max_iter, 20)
            if nrep != self.config.nrep or max_iter != self.config.max_iter:
                self.logger.info(
                    f"MD-Cat adaptive scaling for seq_len={seq_len}: "
                    f"nrep={nrep} (cfg={self.config.nrep}), "
                    f"max_iter={max_iter} (cfg={self.config.max_iter})"
                )

        executable = self._get_executable()
        output_tree = self.work_dir / "mdcat_output.tre"

        cmd = [
            executable,
            "-i",
            str(tree_file),
            "-t",
            str(constraints_file),
            "-l",
            str(seq_len),
            "-k",
            str(self.config.ncat),
            "-p",
            str(nrep),
            "--maxIter",
            str(max_iter),
            "-o",
            str(output_tree),
        ]

        # 时间方向与日期格式
        if self.config.backward_time:
            cmd.append("-b")
        if self.config.as_date:
            cmd.append("-d")

        # 随机种子
        if self.common_config.seed is not None:
            cmd.extend(["--randSeed", str(self.common_config.seed)])

        # 详细输出
        if self.common_config.verbose:
            cmd.append("-v")

        # 置信区间（bootstrap 同样按数据量缩放；≥20000 位点时经验上极慢，默认关闭）
        ci_nboots = self.config.ci_nboots
        if ci_nboots > 0:
            if seq_len >= 20000:
                ci_nboots = 0
                self.logger.info(
                    f"MD-Cat CI disabled for seq_len={seq_len} (bootstrap too slow); "
                    "set ci_nboots manually to override"
                )
            elif seq_len >= 5000:
                ci_nboots = min(ci_nboots, 5)
            if ci_nboots > 0:
                if ci_nboots != self.config.ci_nboots:
                    self.logger.info(
                        f"MD-Cat CI bootstrap scaled for seq_len={seq_len}: "
                        f"ci_nboots={ci_nboots} (cfg={self.config.ci_nboots})"
                    )
                # Upstream `md_cat.py` parses `--CI` as a **single** argparse value
                # that it then splits on whitespace:
                #     tokens = args["CI"].split()
                #     nboots, p_lower, p_upper = tokens[0], tokens[1], tokens[2]
                # Three bare argv tokens would make argparse exit with
                # "unrecognized arguments", and a comma-joined token makes
                # `int()` raise ValueError. The only accepted form is therefore one
                # whitespace-separated string, e.g. ``--CI "100 0.025 0.975"``.
                cmd.extend(
                    [
                        "--CI",
                        f"{ci_nboots} {self.config.ci_plower} {self.config.ci_pupper}",
                    ]
                )
                self.logger.debug(
                    f"MD-Cat CLI CI: nboots={ci_nboots}, p_lower={self.config.ci_plower}, "
                    f"p_upper={self.config.ci_pupper}"
                )

        # 时间参数：与直接 Python 导入路径保持一致，将 root_time/leaf_time
        # 透传给 md-cat CLI（直接导入路径通过 root_time/leaf_time/pseudo 参数传递）。
        # 注意：上游 md_cat.py CLI 并不接受 --pseudo（pseudocount 只存在于 Python
        # API 中，实测 `mdcat -h` 无该标志），因此 subprocess 路径不传该参数。
        if self.config.backward_time:
            tR = 1 if not self.config.as_date else None
            tL = 0 if not self.config.as_date else None
        else:
            tR = 0 if not self.config.as_date else None
            tL = 1 if not self.config.as_date else None
        if tR is not None:
            cmd.extend(["--rootTime", str(tR)])
        if tL is not None:
            cmd.extend(["--leafTime", str(tL)])

        # 注释级别
        cmd.extend(["--annotate", str(self.config.annotate_level)])

        self.logger.info(f"Running MD-Cat via subprocess: {' '.join(cmd)}")
        start_time = time.time()

        # 设置合理的超时：默认 30 分钟，大数据下仍然过久时由调用方决定
        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 1800
        )
        result = runner.run(cmd)

        if result.returncode != 0:
            raise ExecutionError(
                f"MD-Cat execution failed (returncode={result.returncode}): {result.stderr}"
            )

        # md-cat 默认也可能输出到 [input].emDate；若指定的 -o 未生成，尝试兜底
        if not output_tree.exists():
            fallback = Path(str(tree_file) + ".emDate")
            if fallback.exists():
                fallback.rename(output_tree)
                self.logger.info(f"Renamed md-cat fallback output to {output_tree}")

        if not output_tree.exists():
            raise ExecutionError(
                f"MD-Cat did not produce expected output tree: {output_tree}"
            )

        elapsed = time.time() - start_time
        self._input_files["output_tree"] = output_tree
        self._input_files["execution_time"] = elapsed

        self.logger.success(f"MD-Cat execution completed in {elapsed:.2f}s")
        return True

    def parse_results(self) -> DatingResult:
        """解析 MD-Cat 输出"""
        if "output_tree" not in self._input_files:
            # 未执行（或执行失败）就到了 parse：保持原有的 ResultParsingError，
            # 而不是变成 KeyError。
            raise ResultParsingError("MD-Cat output tree not found")
        output_tree = self._require_input_path("output_tree")

        if not output_tree.exists():
            raise ResultParsingError("MD-Cat output tree not found")

        with open(output_tree, "r") as f:
            dated_tree_newick = f.read().strip()

        # 解析节点年龄
        node_ages, parse_warnings = self._parse_tree_ages(dated_tree_newick)

        # 获取执行元数据（int/bool 之外还会追加耗时与对数似然这样的 float）
        metadata: Dict[str, Any] = {
            "ncat": self.config.ncat,
            "replicates": self.config.nrep,
            "backward_time": self.config.backward_time,
            "use_direct_import": self._use_direct_import,
        }

        if "execution_time" in self._input_files:
            metadata["execution_time"] = self._require_input_number("execution_time")
        if "log_likelihood" in self._input_files:
            metadata["log_likelihood"] = self._require_input_number("log_likelihood")

        result = DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=self._require_input_number("execution_time", 0.0),
            metadata=metadata,
        )

        for w in parse_warnings:
            result.add_warning(w)

        self.logger.success(f"Parsed MD-Cat results: {len(node_ages)} node ages")
        return result

    def _resolved_ci_type(self) -> "CIType":
        """C-52: derive interval semantics from the configured quantiles, not hard-coded.

        Only the default (0.025, 0.975) quantiles correspond to a 95% CI; any other
        user-chosen quantiles are labelled RANGE so a 90%/99% interval is not printed
        as _CI95 in the comparison table.
        """
        plower = getattr(self.config, "ci_plower", None)
        pupper = getattr(self.config, "ci_pupper", None)
        try:
            if plower is not None and pupper is not None:
                if (
                    abs(float(plower) - 0.025) < 1e-9
                    and abs(float(pupper) - 0.975) < 1e-9
                ):
                    return CIType.CI95
                return CIType.RANGE
        except (TypeError, ValueError):
            return CIType.RANGE
        return CIType.RANGE

    def _parse_tree_ages(
        self, newick: str
    ) -> Tuple[Dict[str, NodeAgeEstimate], List[Any]]:
        """
        从 MD-Cat 输出树解析年龄，通过 MRCA 拓扑匹配校准点

        MD-Cat 格式：[t=0.5261,mu=0.0764]
        如果启用 --CI: [t=0.5261,mu=0.0764,CI={0.4,0.6}]

        注意：MD-Cat CLI 在 subprocess 模式下通常不在输出树中写入 [t=...] 注解，
        而是输出标准 chronogram（年龄编码在分支长度中）。当输出中不存在注解时，
        直接走 chronogram 距离解码，避免无意义回退并消除冗余警告。

        Returns:
            (node_ages, warnings)
        """
        node_ages: Dict[str, NodeAgeEstimate] = {}
        warnings: List[Any] = []

        # 快速检测：输出是否为不带注解的 chronogram
        if "[t=" not in newick:
            chronogram_ages = self._decode_mdcat_chronogram(newick)
            if chronogram_ages:
                self.logger.info(
                    f"MD-Cat: decoded {len(chronogram_ages)} ages from chronogram branch lengths"
                )
            return chronogram_ages, warnings

        # C-52: numeric class must accept a sign and scientific notation (MD-Cat
        # ages/rates are often printed like 4e-02); the old [\d.]+ silently missed
        # those and degraded whole CI annotations to CIType.NONE.
        _num = r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"
        pattern_basic = rf"\[t=({_num})(?:,mu={_num})?\]"
        pattern_ci = rf"\[t=({_num})(?:,mu={_num})?,CI=\{{({_num}),({_num})\}}\]"
        # C-52: warn when annotations mention CI= but no well-formed block matched,
        # instead of silently dropping the interval.
        if "CI=" in newick and not re.search(r",CI=\{", newick):
            self.logger.warning(
                "MD-Cat: node annotations mention 'CI=' but no well-formed "
                "CI={lower,upper} block was recognized; confidence intervals may be "
                "silently dropped from the parsed ages."
            )

        # 尝试通过 MRCA 拓扑匹配
        matched_by_topology = False
        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)

            for cal in self._calibrations:
                if not cal.mrca_leaf_pair:
                    continue

                tip1, tip2 = cal.mrca_leaf_pair
                try:
                    mrca = tree.get_common_ancestor(tip1, tip2)
                    node_str = mrca.write(format=9)

                    match = re.search(pattern_ci, node_str)
                    if match:
                        node_ages[cal.name] = NodeAgeEstimate(
                            mean_age=float(match.group(1)),
                            ci_lower=float(match.group(2)),
                            ci_upper=float(match.group(3)),
                            ci_type=self._resolved_ci_type(),
                        )
                        matched_by_topology = True
                        continue

                    match = re.search(pattern_basic, node_str)
                    if not match:
                        match = re.search(pattern_basic, mrca.name or "")
                    if match:
                        node_ages[cal.name] = NodeAgeEstimate(
                            mean_age=float(match.group(1)), ci_type=CIType.NONE
                        )
                        matched_by_topology = True
                except Exception as e:
                    self.logger.debug(
                        f"Could not extract age for {cal.name} by MRCA: {e}"
                    )
        except Exception:
            pass

        # 拓扑匹配失败时，不使用索引匹配（索引顺序与校准列表不对应，
        # 会将年龄赋给错误节点），改为逐校准按标签名匹配
        if not matched_by_topology:
            matches_ci = list(re.finditer(pattern_ci, newick))
            matches_basic = list(re.finditer(pattern_basic, newick))
            matches = matches_ci if matches_ci else matches_basic

            if matches and self._calibrations:
                warnings.append(
                    SemanticDegradationWarning(
                        "MD-Cat: MRCA topology matching failed for some calibrations; "
                        f"attempted label-based matching from {len(matches)} annotations."
                    )
                )
                # B-27: warnings[-1] is an exception OBJECT; the logger needs a
                # string, otherwise _sanitize_message does path_pattern.sub(obj) ->
                # TypeError, i.e. the code that reports degradation crashes itself.
                self.logger.warning(str(warnings[-1]))
                for cal in self._calibrations:
                    for match in matches:
                        annotation_context = newick[
                            max(0, match.start() - 80) : match.end() + 20
                        ]
                        if cal.name and cal.name.lower() in annotation_context.lower():
                            age = float(match.group(1))
                            if len(match.groups()) >= 3:
                                node_ages[cal.name] = NodeAgeEstimate(
                                    mean_age=age,
                                    ci_lower=float(match.group(2)),
                                    ci_upper=float(match.group(3)),
                                    ci_type=self._resolved_ci_type(),
                                )
                            else:
                                node_ages[cal.name] = NodeAgeEstimate(
                                    mean_age=age, ci_type=CIType.NONE
                                )
                            break
                    else:
                        self.logger.warning(
                            f"MD-Cat: Could not match calibration '{cal.name}' "
                            f"to any annotation in the output tree. Skipping."
                        )

        # ------------------------------------------------------------------
        # chronogram 回退：MD-Cat 在 use_direct_import=True 模式下输出标准
        # Newick chronogram（年龄编码在分支长度中），并不产生 [t=...,mu=...] 注解。
        # 当注解匹配未能提取到任何节点时，按此回退解码：
        #   - 加载 Newick 为 chronogram（backward_time=True：根=0、向叶递增，
        #     所有叶深度 = 根约束上界）；
        #   - 节点年龄 = max_tip_depth − distance_from_root()（Ma）；
        #   - 键优先取校准点名：内部节点 label 等于某校准点名，或该节点是某
        #     校准 mrca_leaf_pair 的 MRCA；否则回退到 node.label / node_<id>。
        # 这避免 MD-Cat 分析成功却被报告为 "Nodes dated: 0"。
        # ------------------------------------------------------------------
        if not node_ages:
            chronogram_ages = self._decode_mdcat_chronogram(newick)
            if chronogram_ages:
                fallback_msg = (
                    "MD-Cat: annotation-based parsing extracted no ages; "
                    f"fell back to chronogram distance decoding ({len(chronogram_ages)} nodes)."
                )
                warnings.append(SemanticDegradationWarning(fallback_msg))
                self.logger.info(fallback_msg)
                node_ages.update(chronogram_ages)

        return node_ages, warnings

    def _decode_mdcat_chronogram(self, newick: str) -> Dict[str, "NodeAgeEstimate"]:
        """对 MD-Cat 输出的 chronogram 做全树解码（分支长度编码年龄）。

        MD-Cat 在 use_direct_import=True 下输出标准 Newick，年龄编码在分支长度
        中（backward_time=True：根=0、向叶递增）。节点年龄由
        ``max_tip_depth − distance_from_root()`` 得到（单位 Ma）。
        优先以校准点名作键（节点 label 命中校准点名，或节点为某校准的 MRCA）；
        否则回退到 ``node.label`` 与 ``node_<id>``。
        """
        from ..models.results import CIType, NodeAgeEstimate

        node_ages: Dict[str, "NodeAgeEstimate"] = {}
        try:
            from dendropy import Tree as _DendropyTree

            tree = _DendropyTree.get(
                data=newick, schema="newick", preserve_underscores=True
            )
        except Exception as e:
            self.logger.debug(f"MD-Cat chronogram fallback could not load Newick: {e}")
            return node_ages

        try:
            if tree.is_rooted is None or tree.is_rooted is False:
                tree.is_rooted = True
        except Exception:
            pass

        try:
            max_tip_depth = max(
                float(n.distance_from_root()) for n in tree.leaf_node_iter()
            )
        except Exception as e:
            self.logger.debug(
                f"MD-Cat chronogram fallback could not compute tip depth: {e}"
            )
            return node_ages

        # 建立 label(校准点名) → calibration 的映射，便于按节点 label 找回校准点名。
        label_to_calname: Dict[str, str] = {}
        for cal in self._calibrations:
            if cal.name:
                label_to_calname[cal.name.lower()] = cal.name

        # 计算各校准 MRCA 的 node id（用全名在树中匹配——mdcat 输出保留全名）。
        mrca_nodeid_to_calname: Dict[int, str] = {}
        for cal in self._calibrations:
            if not cal.name or not cal.mrca_leaf_pair:
                continue
            taxa = []
            for tip_name in cal.mrca_leaf_pair:
                taxon = tree.taxon_namespace.get_taxon(tip_name)
                if taxon is None:
                    break
                taxa.append(taxon)
            if len(taxa) == len(cal.mrca_leaf_pair):
                try:
                    mrca = tree.mrca(taxa=taxa)
                    if mrca is not None:
                        mrca_nodeid_to_calname[id(mrca)] = cal.name
                except Exception:
                    pass

        for node in tree.preorder_node_iter():
            try:
                age = max_tip_depth - float(node.distance_from_root())
            except Exception:
                continue

            # 优先：节点 label == 校准点名（MD-Cat 会把校准点名写入内部节点 label）
            node_label = node.label or ""
            key = label_to_calname.get((node_label or "").lower(), "")

            # 次优：IntNodeX 映射为统一的 internal_node_X
            if not key and node_label:
                key = self._node_label_map.get(node_label, "")

            # 再次：该节点是某校准的 MRCA
            if not key:
                key = mrca_nodeid_to_calname.get(id(node), "")

            # 过滤叶节点：只保留内部/校准节点
            if not node.child_nodes():
                continue

            # 兜底：node.label → node_<id>
            if not key:
                key = node_label if node_label else f"node_{id(node)}"

            if key and key not in node_ages:
                node_ages[key] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)

        return node_ages


# 注册适配器
DatingMethodRegistry.register("mdcat", MDCatMethod)
