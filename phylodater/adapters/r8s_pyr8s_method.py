"""
R8sPyr8sMethod - r8s/pyr8s 定年适配器

支持原版 C 语言 r8s 和 Python 版 pyr8s
根据诊断报告，pyr8s 仅支持 NPRS 算法
"""

from __future__ import annotations

import math
import re
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, TextIO, Tuple, Type, Union

if TYPE_CHECKING:
    from dendropy import Tree

from ..core import DatingMethod, DatingMethodRegistry
from ..core.exceptions import (
    EnvironmentWarning,
    ExecutionError,
    ResultParsingError,
    SemanticDegradationWarning,
)
from ..core.method_interface import InputFileInfo
from ..infrastructure import ProcessRunner, get_logger
from ..infrastructure.configuration import (
    CommonConfig,
    R8sConfig,
    SoftwarePaths,
    ToolConfig,
)
from ..infrastructure.safe_io import safe_writer
from ..models import (
    CalibrationPoint,
    CIType,
    DatingResult,
    FixedAgeConstraint,
    GammaPriorConstraint,
    MaximumAgeConstraint,
    NodeAgeEstimate,
    PhylogeneticTree,
    SkewNormalConstraint,
    SkewTConstraint,
    SoftBoundsConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
    strip_leading_newick_comments,
)


class R8sPyr8sMethod(DatingMethod[R8sConfig]):
    """
    r8s/pyr8s 定年适配器

    特性：
    - 自动检测后端是原版 r8s 还是 pyr8s
    - pyr8s 仅支持 NPRS，强制降级并发出警告
    - 支持原生 Python API 调用（pyr8s）
    """

    def __init__(
        self,
        config: Union[ToolConfig, R8sConfig],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        # 入口保留双形态（历史上两种调用都存在），但在交给基类之前先归一：
        # 基类把 ``self.config`` 记为 ``R8sConfig``，适配器后续只读子配置。
        # ``isinstance`` 与原先的 ``hasattr(config, "r8s")`` 对这两种入参等价：
        # ``ToolConfig`` 必有 ``r8s``，``R8sConfig`` 必无。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            method_config = config.r8s
            if common is None:
                common = config.common
        else:
            method_config = config
        super().__init__(
            method_config, output_dir, software_paths, common_config=common
        )

        self.logger = get_logger()
        # 同 MD-Cat：这份字典既装文件路径也装执行统计量（seq_len / execution_time）。
        self._input_files: Dict[str, InputFileInfo] = {}
        self._calibrations: List[CalibrationPoint] = []
        # 收集需要写入 DatingResult/runtime_metadata 的警告（与 wlogdate/pathd8 同约定）
        self._warnings: List[Any] = []
        self._is_pyr8s = False
        self._use_api = False

    @property
    def method_name(self) -> str:
        return "r8s"

    def _get_executable(self) -> str:
        """获取 r8s/pyr8s 可执行文件路径，优先使用 software_paths 自定义路径。"""
        path = self.get_software_path("r8s_bin") or self.get_software_path("pyr8s_bin")
        return path if path else self.config.r8s_bin

    def validate_environment(self) -> bool:
        """检测 r8s/pyr8s 是否可用"""
        # 检查是否为 pyr8s
        import importlib.util

        if importlib.util.find_spec("pyr8s") is not None:
            try:
                from pyr8s.core import RateAnalysis  # noqa: F401

                self._is_pyr8s = True
                self._use_api = True
                self.software_version = (
                    f"pyr8s {self._package_version('pyr8s') or 'unknown'} (NPRS only)"
                )
                self.logger.info("Detected pyr8s Python package")

                # 强制锁定为 NPRS
                if self.config.method.upper() != "NPRS":
                    self.logger.warning(
                        f"pyr8s only supports NPRS algorithm. "
                        f"Ignoring '{self.config.method}' setting and using NPRS.",
                        category=EnvironmentWarning,
                    )
                    self.config.method = "NPRS"

                return True
            except ImportError:
                pass

        # 检查是否为原版 r8s
        runner = ProcessRunner()
        executable = self._get_executable()
        available = runner.check_executable(executable)

        if available:
            version = runner.get_version(executable, "-v")
            self.software_version = version or f"r8s binary {executable}"
            self.logger.info(f"Detected r8s: {version or 'unknown version'}")
        else:
            self.logger.warning(f"r8s/pyr8s not found: {executable}")

        return available

    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
        生成 r8s/pyr8s 输入文件

        生成 NEXUS 格式文件，包含 TAXA、TREES 和 RATES/R8S 块。
        通过显式 MRCA 定义 + FIXAGE/CONSTRAIN 命令把校准约束真正传递给后端：
          - pyr8s 使用 BEGIN RATES; ... BLFORMAT ... MRCA ... FIXAGE ... 语法
          - 原版 r8s 使用 BEGIN R8S; ... MRCA ... FIXAGE ... 语法
        """
        self._calibrations = calibrations
        #: 每次运行重新收集结果级警告，避免复用实例时上一次的降级记录串进本次结果
        self._warnings = []
        self._warn_pyr8s_advisory_ranges(calibrations)

        # 序列长度
        if alignment_path and alignment_path.exists():
            from ..infrastructure import AlignmentMetadataExtractor

            extractor = AlignmentMetadataExtractor()
            metadata = extractor.extract(alignment_path)
            seq_len = metadata.sequence_length
        else:
            # C-42：无比对时不再凭空使用 1000 而不留痕。该值会写进
            # blformat nsites= / BLFORMAT NSITES= 参与速率与方差换算，必须显式披露，
            # 否则"跑了比对"与"没跑比对"的两种情形会产出看似同源的结果。
            seq_len = 1000
            self.logger.warning(
                "r8s: no alignment provided; assuming sequence length = 1000. "
                "This value is written into blformat nsites= and scales rate/variance "
                "estimates. Provide an alignment for reproducible results.",
                category=SemanticDegradationWarning,
            )

        # 建立安全短名映射，使 NEXUS 中的 tip 名与 MRCA 命令中的名一致。
        # r8s/pyr8s 对长名称/特殊字符敏感，统一用 _sanitize_name 处理。
        tip_names = tree.tip_names
        name_map: Dict[str, str] = {}
        used = set()
        for i, name in enumerate(tip_names):
            base = self._sanitize_name(name)
            safe = base
            suffix = 0
            while safe in used or not safe:
                suffix += 1
                safe = f"{base[:25]}_{suffix}"
            name_map[name] = safe
            used.add(safe)
        self._tip_name_map = name_map

        # 用短名重写树并去掉内部节点标签
        renamed_tree = tree.with_renamed_leaves(name_map).without_internal_labels()

        # 生成 NEXUS 文件
        nexus_path = self.work_dir / "r8s_input.nex"
        block_name = "RATES" if self._is_pyr8s else "R8S"

        with safe_writer(nexus_path) as f:
            f.write("#NEXUS\n")
            f.write("BEGIN TAXA;\n")
            f.write(f"    DIMENSIONS NTAX={len(tip_names)};\n")
            f.write("    TAXLABELS\n")
            for name in tip_names:
                f.write(f"        {name_map[name]}\n")
            f.write("    ;\n")
            f.write("END;\n\n")

            f.write("BEGIN TREES;\n")
            f.write(f"    TREE tree1 = {renamed_tree.newick}\n")
            f.write("END;\n\n")

            f.write(f"BEGIN {block_name};\n")

            # 分支长度格式设置
            # pyr8s 沿用其大写 BLFORMAT 语法（已核验通路，保持不动）。
            # 原版 r8s 依 Sanderson 1.7 手册：blformat 必须是 r8s 块内首条命令
            #   blformat lengths=total|persite nsites=nnnn ultrametric=no|yes;
            # 旧的 blength/method=/algorithm=/smoothing= 四条命令上游根本不存在，
            # 会让方法/算法/平滑三项配置全部送不进 r8s（A-6）。
            if self._is_pyr8s:
                f.write(f"    BLFORMAT LENGTHS=PERSITE NSITES={seq_len};\n")
            else:
                f.write(
                    f"    blformat lengths=persite nsites={seq_len} ultrametric=no;\n"
                )
                # 手册：set smoothing=<real>;（仅当用户显式给出时对 PL/NPRS 生效）
                if self.config.smoothing is not None:
                    f.write(f"    set smoothing={self.config.smoothing};\n")

            # 根节点校准：用跨根的两个叶子定义 MRCA
            root_cal = next((c for c in calibrations if c.is_root_node), None)
            if root_cal and root_cal.age_constraint is None:
                self.logger.warning(
                    f"Skipping root calibration '{root_cal.name}' in r8s: "
                    "no age_constraint（无年龄界可写）"
                )
                root_cal = None
            if root_cal:
                rt1, rt2 = self._get_root_representative_leaves(tree, name_map)
                f.write(f"    MRCA Root {rt1} {rt2};\n")
                self._write_r8s_constraint(f, root_cal, "Root")
                self._warn_constraint_degradation(root_cal)

            # 内部校准点
            for cal in calibrations:
                if cal.is_root_node:
                    continue
                # B-20：mrca_leaf_pair 是 Optional，根校准/缺失叶对都合法，
                # 必须在下标前显式守卫，否则对 None 下标会 TypeError。
                pair = cal.mrca_leaf_pair
                if not pair or len(pair) < 2:
                    self.logger.warning(
                        f"Skipping calibration '{cal.name}' in r8s: "
                        f"missing or short mrca_leaf_pair ({pair!r})"
                    )
                    continue
                if cal.age_constraint is None:
                    # 与上面“缺叶对就跳过”同一条口径：没有年龄界的校准点写不出
                    # FIXAGE/CONSTRAIN，而 MRCA 行一旦写出去就会被 r8s 当成一个
                    # 有名无实的约束节点，所以必须在写 MRCA 之前就跳过。
                    self.logger.warning(
                        f"Skipping calibration '{cal.name}' in r8s: no age_constraint"
                    )
                    continue

                t1 = name_map.get(pair[0], self._sanitize_name(pair[0]))
                t2 = name_map.get(pair[-1], self._sanitize_name(pair[-1]))
                label = self._sanitize_name(cal.name)
                f.write(f"    MRCA {label} {t1} {t2};\n")
                self._write_r8s_constraint(f, cal, label)
                self._warn_constraint_degradation(cal)

            # 定年与输出
            if self._is_pyr8s:
                f.write("    DIVTIME METHOD=NPRS ALGORITHM=POWELL;\n")
                f.write("    DESCRIBE PLOT=CHRONOGRAM;\n")
            else:
                # 手册：divtime method=LF|NPRS|PL algorithm=POWELL|TN|QNEWT;
                f.write(
                    f"    divtime method={self.config.method} "
                    f"algorithm={self.config.algorithm};\n"
                )
                # 手册：showage shownamed=yes 只列出被 MRCA 命令命名的节点年龄
                f.write("    showage shownamed=yes;\n")
                # 手册：chrono_description 才打印可导入的 Nexus 树描述（分支长度=时间）；
                # plot=chronogram 只输出 ASCII 字符图，parse 端拿不到可解析的定年树（A-6）。
                f.write("    describe plot=chrono_description;\n")
            f.write("END;\n")

        self._input_files = {"nexus": nexus_path, "seq_len": seq_len}
        self.logger.info(f"Generated r8s/pyr8s NEXUS file: {nexus_path}")

        return self._input_files

    def _get_root_representative_leaves(
        self, tree: PhylogeneticTree, name_map: Dict[str, str]
    ) -> Tuple[str, str]:
        """从根节点的两个子树各取一片叶子"""
        try:
            from ete3 import Tree

            et = Tree(strip_leading_newick_comments(tree.newick), format=1)
            root = et.get_tree_root()
            children = [ch for ch in root.children if not ch.is_leaf()]
            if len(children) < 2:
                children = root.children
            leaves = []
            for ch in children[:2]:
                leaf = next(ch.iter_leaves(), None)
                if leaf is not None:
                    leaves.append(
                        name_map.get(leaf.name, self._sanitize_name(leaf.name))
                    )
            if len(leaves) == 2:
                return leaves[0], leaves[1]
        except Exception as e:
            self.logger.warning(f"Could not pick root representative leaves: {e}")

        # 回退：使用前两个 tip
        names = list(name_map.values())[:2]
        return names[0], names[1]

    def _write_r8s_constraint(
        self, f: TextIO, cal: CalibrationPoint, label: str
    ) -> None:
        """将单个校准约束写成 MRCA 后的 FIXAGE/CONSTRAIN 命令。

        两个调用点都已在写 MRCA 行前拦掉 ``age_constraint is None``，走到这里
        还是 None 就是调用顺序被破坏了，必须报错而不是静默补 0。
        """
        constraint = cal.age_constraint
        if constraint is None:
            raise ExecutionError(
                f"r8s/pyr8s: calibration '{cal.name}' has no age_constraint, so no "
                "FIXAGE/CONSTRAIN command can be written for it"
            )
        target = "pyr8s" if self._is_pyr8s else "r8s"
        formatted = constraint.to_software_format(target, label=label)
        for line in formatted.splitlines():
            line = line.strip()
            if not line:
                continue
            if self._is_pyr8s:
                line = self._pyr8s_integer_ages(line)
            f.write(f"    {line}\n")

    # C-40：pyr8s 只接受整数年龄。这里显式半进位取整（避免 Python round 的
    # 银行家舍入），兼容 .5 这类省略整数位的写法，并在确实改变取值时发出
    # 降级警告而不是静默改数；(?<![\w.]) 保证不改 Node1.2 这类含点的节点名。
    _PYR8S_AGE_NUM = re.compile(r"(?<![\w.])\d*\.\d+(?![\w.])")

    def _pyr8s_integer_ages(self, line: str) -> str:
        """pyr8s 的 AGE/MINAGE/MAXAGE 只接受整数，把浮点年龄显式取整并留痕。"""

        def _convert(match: "re.Match[str]") -> str:
            text = match.group(0)
            value = float(text)
            rounded = math.floor(value + 0.5)
            if abs(value - rounded) > 1e-12:
                self.logger.warning(
                    f"pyr8s requires integer ages: {text} -> {rounded} "
                    f"(precision lost; AGE/MINAGE/MAXAGE accept integers only).",
                    category=SemanticDegradationWarning,
                )
            return str(rounded)

        return self._PYR8S_AGE_NUM.sub(_convert, line)

    def _warn_constraint_degradation(self, cal: CalibrationPoint) -> None:
        """为降级的约束类型发出适配器级警告"""
        constraint = cal.age_constraint
        if isinstance(constraint, SoftLowerBoundConstraint):
            self.logger.warning(
                f"Soft lower bound degraded to hard min_age for r8s: {cal.name}",
                category=SemanticDegradationWarning,
            )
        elif isinstance(constraint, MaximumAgeConstraint):
            self.logger.warning(
                f"Maximum age degraded to hard max_age for r8s: {cal.name}",
                category=SemanticDegradationWarning,
            )
        elif isinstance(
            constraint,
            (
                SoftBoundsConstraint,
                GammaPriorConstraint,
                SkewNormalConstraint,
                SkewTConstraint,
            ),
        ):
            self.logger.warning(
                f"{type(constraint).__name__} degraded to uniform range for r8s: {cal.name}",
                category=SemanticDegradationWarning,
            )

    def _warn_pyr8s_advisory_ranges(self, calibrations: List[CalibrationPoint]) -> None:
        """pyr8s 后端下，非点校准（CONSTRAIN 区间）不会被 NPRS 强制满足——必须预先披露。

        已用随环境提供的 pyr8s 0.3.1 直接实验确认（见
        ``scripts/probe_pyr8s_constraints.py``）：对内部节点写
        ``constrain taxon=X min_age=68; max_age=82;`` 时 NPRS 仍可以返回 54.6 Ma；
        而同一个节点改用 ``fixage`` 时年龄被精确固定。也就是说 pyr8s 的 CONSTRAIN
        是搜索中的势垒惩罚（会被逐轮松弛），不是硬约束。

        本项目的纪律是"绝不对用户静默丢约束"，因此在写输入文件时就把这件事说出来，
        并在 ``parse_results()`` 后逐个区间实测回报是否被满足。
        """
        if not self._is_pyr8s:
            return
        advisory = [
            cal.name
            for cal in calibrations
            if not isinstance(cal.age_constraint, FixedAgeConstraint)
        ]
        if not advisory:
            return
        msg = (
            "pyr8s backend: CONSTRAIN ranges are ADVISORY, not enforced by NPRS "
            f"(measured: pyr8s 0.3.1 returns ages outside the requested bounds for "
            f"{', '.join(advisory)}). Only FIXAGE (point) calibrations are honoured. "
            "Use mcmctree/treepl if the requested range must be respected, or convert "
            "these calibrations to fixed ages. The post-run check below reports whether "
            "each range was actually satisfied."
        )
        self.logger.warning(msg, category=SemanticDegradationWarning)
        self._warnings.append(SemanticDegradationWarning(msg))

    def _verify_requested_constraints(
        self, node_ages: Dict[str, NodeAgeEstimate]
    ) -> None:
        """把返回的节点年龄与用户请求的约束逐条对账，未满足则写进结果级警告。

        只对能写成区间的约束类型（fixed / uniform / soft_bounds / maximum /
        soft_lower）检查；概率型约束已在渲染时被降级，不在此重复评判。
        """
        for cal in self._calibrations:
            constraint = cal.age_constraint
            estimate = node_ages.get(cal.name)
            if constraint is None or estimate is None:
                continue
            low: Optional[float] = None
            high: Optional[float] = None
            if isinstance(constraint, FixedAgeConstraint):
                low = high = constraint.fixed_age
            elif isinstance(constraint, (UniformAgeConstraint, SoftBoundsConstraint)):
                low, high = constraint.min_age, constraint.max_age
            elif isinstance(constraint, MaximumAgeConstraint):
                high = constraint.max_age
            elif isinstance(constraint, SoftLowerBoundConstraint):
                low = constraint.min_age
            else:
                continue

            age = float(estimate.mean_age)
            # pyr8s 只接受整数年龄，因此允许 1 Ma 的取整边界
            tolerance = 1.0 + 1e-6
            broken: List[str] = []
            if low is not None and age < low - tolerance:
                broken.append(f"younger than the requested minimum {low} Ma")
            if high is not None and age > high + tolerance:
                broken.append(f"older than the requested maximum {high} Ma")
            if broken:
                msg = (
                    f"{self.method_name}: calibration '{cal.name}' was NOT satisfied "
                    f"by the backend: node age {age:.2f} Ma "
                    f"({' and '.join(broken)}). Treat this estimate as inconsistent "
                    "with the requested calibration."
                )
                self.logger.warning(msg, category=SemanticDegradationWarning)
                self._note_result_warning(msg)

    def _sanitize_name(self, name: str) -> str:
        """安全化名称"""
        safe = re.sub(r"[^A-Za-z0-9_]", "_", name)
        return safe[:30]

    def execute(self) -> bool:
        """执行 r8s/pyr8s 分析"""
        if self._is_pyr8s and self._use_api:
            return self._execute_with_api()
        return self._execute_with_subprocess()

    def _execute_with_api(self) -> bool:
        """使用 pyr8s Python API 执行

        通过 pyr8s.parse.from_file_nexus 解析 RATES 块并运行 NPRS，
        这样校准约束会被真正读取，而不是被忽略。
        """
        try:
            from pyr8s.parse import from_file_nexus

            start_time = time.time()
            analysis = from_file_nexus(str(self._require_input_path("nexus")), run=True)
            self._input_files["execution_time"] = time.time() - start_time
            chronogram = analysis.results.chronogram

            output_tree_path = self.work_dir / "r8s_output.tre"
            chronogram.write_to_path(
                str(output_tree_path),
                schema="newick",
                suppress_internal_node_labels=True,
            )

            self._input_files["output_tree"] = output_tree_path
            self.logger.success("pyr8s API execution completed")
            return True

        except Exception as e:
            raise ExecutionError(f"pyr8s API execution failed: {e}")

    def _execute_with_subprocess(self) -> bool:
        """使用子进程执行（原版 r8s）"""
        nexus_file = self._require_input_path("nexus")

        # A-6：手册要求批处理模式加 -b，并用 -f 指定输入文件；否则 r8s
        # "会耐心地等待键盘输入"（终端下挂到超时；管道下可能一条命令不跑返回 0）。
        cmd = [self._get_executable(), "-f", str(nexus_file), "-b"]

        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 3600
        )
        start_time = time.time()
        result = runner.run(cmd)
        self._input_files["execution_time"] = time.time() - start_time

        if result.returncode != 0:
            raise ExecutionError(
                f"r8s execution failed (returncode={result.returncode}): {result.stderr}"
            )

        # 保存输出
        output_path = self.work_dir / "r8s_output.txt"
        with safe_writer(output_path) as f:
            f.write(result.stdout)

        self._input_files["output"] = output_path
        self.logger.success("r8s execution completed")
        return True

    def parse_results(self) -> DatingResult:
        """解析 r8s/pyr8s 输出"""
        if self._is_pyr8s and self._use_api:
            return self._parse_api_results()
        return self._parse_subprocess_results()

    def _parse_api_results(self) -> DatingResult:
        """解析 pyr8s API 执行结果"""
        if "output_tree" not in self._input_files:
            # 没执行过就到 parse：保持原有的 ResultParsingError，而不是 KeyError。
            raise ResultParsingError("pyr8s output tree not found")
        output_tree_path = self._require_input_path("output_tree")

        if not output_tree_path.exists():
            raise ResultParsingError("pyr8s output tree not found")

        with open(output_tree_path, "r") as f:
            dated_tree_newick = f.read().strip()

        # 将短名叶节点恢复为原始全名
        dated_tree_newick = self._restore_leaf_names(dated_tree_newick)

        try:
            from dendropy import Tree

            chronogram = Tree.get(
                data=dated_tree_newick, schema="newick", preserve_underscores=True
            )
        except Exception as e:
            raise ResultParsingError(f"Failed to parse pyr8s chronogram: {e}")

        node_ages = self._extract_ages_from_chronogram(chronogram)
        self._validate_node_ages(node_ages)
        self._verify_requested_constraints(node_ages)

        self.logger.success(f"Parsed pyr8s API results: {len(node_ages)} node ages")

        return DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=self._require_input_number("execution_time", 0.0),
            # C-41：被排除的节点等降级信息必须随结果落盘，而不是只留在日志里
            warnings=list(self._warnings),
        )

    @staticmethod
    def _describe_tree_node(node: Any) -> str:
        """给一个 dendropy 内部节点一个真能认出"是哪个节点"的名字。

        同时给出**节点标签**（如果有）与**后代叶集合**：C-41 要求点名，而
        ``internal_node_<遍历序号>`` 那种位置键、或只有一个上游残留标签，都不足以让
        人在比较表里找回这一行。叶太多时截断并给出总数，避免刷屏。
        """
        try:
            label = str(getattr(node, "label", None) or "").strip()
        except Exception:
            label = ""
        try:
            leaves = sorted(
                str(leaf.taxon.label)
                for leaf in node.leaf_iter()
                if leaf.taxon is not None
            )
        except Exception:
            leaves = []
        shown = ", ".join(leaves[:5]) + (" …" if len(leaves) > 5 else "")
        parts = [f"label={label!r}"] if label else []
        parts.append(f"leaves=[{shown}] ({len(leaves)} leaves)")
        return ", ".join(parts)

    def _note_result_warning(
        self,
        msg: str,
        category: Type[Warning] = SemanticDegradationWarning,
    ) -> Warning:
        """响亮地报出一条警告，并把它挂到本次运行的结果上（随产物落盘）。

        只在日志里响一声是不够的：日志会滚走，而"少了哪些节点"必须留在产物里
        （``pipeline`` 把 ``result.warnings`` 写进 summary/runtime_metadata）。
        """
        self.logger.warning(msg, category=category)
        warning = category(msg)
        self._warnings.append(warning)
        return warning

    def _extract_ages_from_chronogram(
        self, chronogram: "Tree"
    ) -> Dict[str, NodeAgeEstimate]:
        """从 pyr8s chronogram 中按 MRCA 拓扑提取节点年龄

        chronogram 的分支长度表示绝对时间跨度。一个节点的年龄等于该节点
        到其所有后代叶子路径长度的最大值（即到最近叶子的剩余时间）。
        根节点年龄即树高（到最远叶子的距离）。

        C-41 纪律：本函数**不静默丢节点**。年龄算成 ``<= 0`` 的未命名内部节点仍然不
        会被写成一行年龄（非正年龄不是可用的估计），但会被点名记进
        ``self._warnings`` 并随 ``DatingResult.warnings`` 落盘——"比较表少了行"这件事
        必须能在产物里看见。
        """
        node_ages: Dict[str, NodeAgeEstimate] = {}

        # 确保 dendropy 按有根树处理，否则 mrca() 会错误地返回 seed_node
        try:
            if chronogram.is_rooted is None or chronogram.is_rooted is False:
                chronogram.is_rooted = True
        except Exception:
            pass

        def _max_descendant_distance(node: Any) -> float:
            """节点到其最远后代叶子的距离"""
            try:
                node_root_dist = float(node.distance_from_root())
                return max(
                    float(leaf.distance_from_root()) - node_root_dist
                    for leaf in node.leaf_iter()
                )
            except Exception:
                return 0.0

        for cal in self._calibrations:
            if not cal.name:
                continue

            if cal.is_root_node:
                node_ages[cal.name] = NodeAgeEstimate(
                    mean_age=_max_descendant_distance(chronogram.seed_node),
                    ci_type=CIType.NONE,
                )
                continue

            pair = cal.mrca_leaf_pair
            if not pair or len(pair) < 2:
                continue

            # 使用与输入树一致的名称定位 MRCA。parse_results 中已先调用
            # _restore_leaf_names 把输出树恢复为原始全名，因此优先用全名
            # 查找；若找不到再回退到短名（兼容未改名或中间流程未恢复的情况）。
            taxa = []
            for tip in pair:
                taxon = chronogram.taxon_namespace.get_taxon(tip)
                if taxon is None and tip in self._tip_name_map:
                    taxon = chronogram.taxon_namespace.get_taxon(
                        self._tip_name_map[tip]
                    )
                if taxon is None:
                    taxon = chronogram.taxon_namespace.get_taxon(
                        self._sanitize_name(tip)
                    )
                if taxon is None:
                    break
                taxa.append(taxon)
            if len(taxa) != len(pair):
                continue

            try:
                mrca = chronogram.mrca(taxa=taxa)
                age = _max_descendant_distance(mrca)
            except Exception:
                continue

            node_ages[cal.name] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)

        # 补充未命名内部节点，使 comparison_table 与其他方法一致
        try:
            calibrated_node_ids = set()
            for cal in self._calibrations:
                if not cal.name:
                    continue
                if cal.is_root_node:
                    calibrated_node_ids.add(id(chronogram.seed_node))
                    continue
                pair = cal.mrca_leaf_pair
                if not pair or len(pair) < 2:
                    continue
                taxa = []
                for tip in pair:
                    taxon = chronogram.taxon_namespace.get_taxon(tip)
                    if taxon is None and tip in getattr(self, "_tip_name_map", {}):
                        taxon = chronogram.taxon_namespace.get_taxon(
                            self._tip_name_map[tip]
                        )
                    if taxon is None:
                        taxon = chronogram.taxon_namespace.get_taxon(
                            self._sanitize_name(tip)
                        )
                    if taxon is None:
                        break
                    taxa.append(taxon)
                if len(taxa) != len(pair):
                    continue
                try:
                    mrca = chronogram.mrca(taxa=taxa)
                    if mrca is not None:
                        calibrated_node_ids.add(id(mrca))
                except Exception:
                    pass

            internal_idx = 1
            dropped_nodes: List[str] = []
            for node in chronogram.preorder_node_iter():
                if node.is_leaf():
                    continue
                if id(node) in calibrated_node_ids:
                    continue
                age = _max_descendant_distance(node)
                if age <= 0:
                    # C-41：旧写法是 ``continue``——节点就这样从结果里消失了，不留痕迹
                    # （比较表少一行、日志报"分析成功"，没人知道少了谁）。非正年龄本身
                    # 不该被当成年龄印进表里，所以仍然不产出该行，但必须**点名**披露：
                    # 哪个节点（按标签/后代叶集合识别，而不是遍历序号）、算出了多少。
                    description = self._describe_tree_node(node)
                    dropped_nodes.append(f"{description} -> age={age:g}")
                    continue
                while f"internal_node_{internal_idx}" in node_ages:
                    internal_idx += 1
                node_ages[f"internal_node_{internal_idx}"] = NodeAgeEstimate(
                    mean_age=age, ci_type=CIType.NONE
                )
                internal_idx += 1

            if dropped_nodes:
                self._note_result_warning(
                    "r8s: "
                    f"{len(dropped_nodes)} unnamed internal node(s) were excluded from "
                    "the results because their computed age was <= 0 (a non-positive "
                    "age is not a usable estimate, but silently dropping nodes would "
                    "make comparison_table.tsv look complete while it is not): "
                    + "; ".join(dropped_nodes)
                    + ". Check the input tree's zero/negative branch lengths or the "
                    "r8s smoothing solution."
                )
        except Exception as e:
            # 整段提取失败 = 所有未命名内部节点都没有年龄。只记 debug 等于什么都不说
            # （与 C-41 同一形态的静默），因此升级为结果级警告。
            self._note_result_warning(
                f"r8s: could not extract unnamed internal node ages "
                f"({e}); comparison_table.tsv will contain calibrated nodes only."
            )

        return node_ages

    def _parse_subprocess_results(self) -> DatingResult:
        """解析原版 r8s 子进程输出（showage shownamed=yes + chrono_description）。

        年龄取自 ``describe plot=chrono_description`` 输出的定年树：按拓扑（校准点的
        mrca_leaf_pair / 根）用与已核验的 pyr8s API 通路完全相同的
        `_extract_ages_from_chronogram` 提取节点年龄，避免旧实现"子串匹配 + 行内第一个
        数字"把 80 Ma 误读成 1.0 Ma（B-25）。取不到可解析且含末端节点的树、或解析不到
        任何年龄时显式 raise，而不是以 PASS 返回空结果。
        """
        if "output" not in self._input_files:
            raise ResultParsingError("r8s output not found")
        output_path = self._require_input_path("output")

        if not output_path.exists():
            raise ResultParsingError("r8s output not found")

        with open(output_path, "r") as f:
            content = f.read()

        newick = self._extract_chrono_description(content)
        if not newick:
            raise ResultParsingError(
                "r8s: could not locate a `chrono_description` tree in the output. "
                "The r8s block must emit `describe plot=chrono_description;` (a Nexus "
                "tree description); `plot=chronogram` only prints an ASCII character "
                "plot. No importable dated tree was produced, so the run is failed."
            )

        dated_tree_newick = self._restore_leaf_names(newick)

        try:
            from dendropy import Tree

            chronogram = Tree.get(
                data=dated_tree_newick,
                schema="newick",
                preserve_underscores=True,
            )
        except Exception as e:
            raise ResultParsingError(
                f"r8s: failed to parse chrono_description tree: {e}"
            )

        if chronogram is None or not list(chronogram.leaf_node_iter()):
            raise ResultParsingError(
                "r8s: the chrono_description tree parsed but contains no terminal nodes."
            )

        node_ages = self._extract_ages_from_chronogram(chronogram)

        # 回退：若拓扑提取为空（例如输出树缺内部信息），改按
        # showage shownamed=yes 的"整字段名 -> 年龄"列对齐解析（不做子串匹配）。
        if not node_ages:
            node_ages = self._parse_showage_named_ages(content)

        self._validate_node_ages(node_ages)

        result = DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=self._require_input_number("execution_time", 0.0),
            # C-41：被排除的节点等降级信息必须随结果落盘，而不是只留在日志里
            warnings=list(self._warnings),
        )

        self.logger.success(f"Parsed r8s results: {len(node_ages)} node ages")
        return result

    def _extract_chrono_description(self, content: str) -> str:
        """从 r8s stdout 取出 ``describe plot=chrono_description`` 打印的定年树 Newick。

        r8s 以 Nexus 树描述形式打印，典型为 ``tree "chronogram" = <newick>;``，可能跨
        行、并带 ``[&R]`` 有根标记。这里定位 ``tree … =``，做括号配平截取到结束分号。
        取不到返回空串（由调用方 raise）。
        """
        for m in re.finditer(r"tree\s+[^\s=]*\s*=", content, re.IGNORECASE):
            start = content.find("(", m.end())
            if start == -1:
                continue
            end = self._match_balanced_paren(content, start)
            if end is None:
                continue
            semi = content.find(";", end + 1)
            # 允许结束括号与分号之间存在根标签/长度（如 ")Root:0.0"），但不宜过长
            if semi == -1 or (semi - end) > 128:
                continue
            newick = content[start : semi + 1].strip()
            if newick.count("(") >= 1 and "," in newick:
                return newick
        return ""

    @staticmethod
    def _match_balanced_paren(text: str, start: int) -> Optional[int]:
        """返回与 text[start]=='(' 匹配的右括号下标，找不到返回 None。"""
        if start >= len(text) or text[start] != "(":
            return None
        depth = 0
        for i in range(start, len(text)):
            c = text[i]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    return i
        return None

    # showage 表里年龄列以普通十进制打印，速率列多为科学计数法（如 1.2e-02）、
    # Type 列是整数——用"纯十进制"这一特征稳定地只挑出年龄列。
    _PLAIN_AGE = re.compile(r"^-?\d+\.\d+$")

    def _parse_showage_named_ages(self, content: str) -> Dict[str, NodeAgeEstimate]:
        """按"整字段名 -> 年龄"列对齐解析 r8s `showage shownamed=yes` 表（B-25 回退）。

        仅当某行第一个字段与校准点名完全相等（不做子串匹配）时命中；年龄取该行中唯一
        的纯十进制数值列。同一名字命中多行且年龄不一致时 raise，绝不"最后写入者胜出"。
        """
        node_ages: Dict[str, NodeAgeEstimate] = {}
        cal_names = {c.name for c in self._calibrations if c.name}
        if not cal_names:
            return node_ages

        for raw in content.splitlines():
            line = raw.strip()
            if not line:
                continue
            fields = line.split()
            name = fields[0]
            if name not in cal_names:
                continue
            age_values = sorted(
                {round(float(t), 9) for t in fields[1:] if self._PLAIN_AGE.match(t)}
            )
            if not age_values:
                continue
            if len(age_values) > 1:
                raise ResultParsingError(
                    f"r8s showage: calibration '{name}' row matched multiple "
                    f"conflicting age columns {age_values}; refusing to guess."
                )
            value = age_values[0]
            # 同一名字命中多行且年龄不一致：raise，绝不"最后写入者胜出"（B-25）。
            if name in node_ages and abs(node_ages[name].mean_age - value) > 1e-9:
                raise ResultParsingError(
                    f"r8s showage: calibration '{name}' matched conflicting ages "
                    f"{node_ages[name].mean_age} and {value} across rows; the r8s "
                    f"output format may have changed. Refusing to guess."
                )
            node_ages[name] = NodeAgeEstimate(mean_age=value, ci_type=CIType.NONE)
        return node_ages

    def _validate_node_ages(self, node_ages: Dict[str, NodeAgeEstimate]) -> None:
        """校验提取到的节点年龄：空结果与全 0 结果都视为解析失败（B-25）。

        r8s/pyr8s 在约束类型不兼容（如 uniform/maximum）时可能成功退出但输出全 0，
        或因输出格式不符而一个年龄都没解析出来。旧实现在 `if not node_ages: return`
        处短路，让"什么都没解析出来"以 PASS 返回空结果；此处两者一律 raise。
        """
        if not node_ages:
            raise ResultParsingError(
                "r8s/pyr8s: no node ages could be parsed from the output. The dated "
                "tree may be missing, empty, or not match the requested calibrations. "
                "Failing loudly instead of reporting success with zero node ages."
            )

        ages = [estimate.mean_age for estimate in node_ages.values()]
        if all(age == 0.0 for age in ages):
            raise ResultParsingError(
                "r8s/pyr8s returned all-zero ages. This usually means the selected "
                "constraint types are not supported by the backend (pyr8s NPRS only "
                "reliably handles fixed ages). Consider using fixed constraints or "
                "the original r8s binary."
            )

    def _restore_leaf_names(self, newick: str) -> str:
        """将输出树中的短名叶节点替换回原始全名，并校正内部校准节点标签。

        prepare_inputs 阶段使用 `_tip_name_map` 把原始全名映射为 ≤30 字符的
        安全短名；本方法建立反向映射，通过 ETE3 解析树并替换叶节点标签。
        此外，pyr8s 可能把内部校准名中的 '-' 替换为 '_'，这里用校准点原始名称
        重新标注 MRCA 节点，保证输出标签与输入配置一致。
        pyr8s API 输出可能以 `[&R] ` 开头，ETE3 format=1 无法直接解析该前缀，
        因此先临时移除、处理完成后再加回。
        """
        tip_name_map = getattr(self, "_tip_name_map", None)
        if not tip_name_map:
            return newick

        reverse_mapping = {short: full for full, short in tip_name_map.items()}

        # pyr8s 输出可能以 [&R] 标记有根树，ETE3 format=1 不支持，先临时移除
        rooted_prefix = "[\u0026R] "
        prefix = ""
        tree_text = newick
        if tree_text.startswith(rooted_prefix):
            prefix = rooted_prefix
            tree_text = tree_text[len(rooted_prefix) :]

        try:
            from ete3 import Tree

            tree = Tree(tree_text, format=1)

            # 1. 恢复叶节点名
            for leaf in tree.get_leaves():
                full_name = reverse_mapping.get(leaf.name)
                if full_name:
                    leaf.name = full_name

            # 2. 用原始校准名重新标注内部 MRCA 节点
            for cal in getattr(self, "_calibrations", []):
                if not cal.name:
                    continue
                if cal.is_root_node:
                    tree.name = cal.name
                    continue
                if not cal.mrca_leaf_pair or len(cal.mrca_leaf_pair) < 2:
                    continue
                try:
                    ancestor = tree.get_common_ancestor(*cal.mrca_leaf_pair)
                    ancestor.name = cal.name
                except Exception:
                    pass

            restored: str = tree.write(format=1)
            return prefix + restored if prefix else restored
        except ImportError:
            self.logger.warning("ETE3 not available, cannot restore r8s leaf names")
            return newick
        except Exception as e:
            self.logger.warning(f"Failed to restore r8s leaf names: {e}")
            return newick


# 注册适配器
DatingMethodRegistry.register("r8s", R8sPyr8sMethod)
