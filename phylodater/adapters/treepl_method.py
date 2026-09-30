"""
TreePLMethod - treePL 定年适配器

实现符合最佳实践的三阶段定年流程：
1. Prime步骤 - 获取最佳优化参数(opt/optad)
2. CV分析 - 获取最佳平滑参数
3. 最终定年 - 使用优化参数和平滑参数进行定年

参考: treePL经验教程.pdf
"""

import re
import secrets
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, TextIO, Tuple, Union

from ..core import DatingMethod, DatingMethodRegistry
from ..core.exceptions import (
    CalibrationError,
    ConvergenceWarning,
    ExecutionError,
    ResultParsingError,
    SemanticDegradationWarning,
)
from ..infrastructure import AlignmentMetadataExtractor, ProcessRunner, get_logger
from ..infrastructure.configuration import (
    CommonConfig,
    SoftwarePaths,
    ToolConfig,
    TreePLConfig,
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


@dataclass
class PrimeResult:
    """Prime步骤结果"""

    opt: int
    optad: int
    success: bool
    message: str = ""
    stdout: str = ""  # 保留原始输出以便错误诊断


@dataclass
class CVResult:
    """CV分析结果"""

    smoothing: float
    chisq: float
    success: bool
    message: str = ""


class TreePLMethod(DatingMethod[TreePLConfig]):
    """
    treePL 定年适配器

    特性：
    - 三阶段流程：Prime → CV → 最终定年
    - 支持randomcv（推荐）和cv两种交叉验证方法
    - 多叶节点 MRCA 支持
    - CV重复运行以检查平滑值稳定性
    - Prime重复运行选择最低opt/optad
    - 自动检测和处理tiny branch
    - 自动调整opt/optad消除警告
    - 自适应CV参数调整
    - 随机种子确保可重复性
    """

    def __init__(
        self,
        config: Union[ToolConfig, TreePLConfig],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        # 入口保留双形态（历史上两种调用都存在），但在交给基类之前先归一：
        # 基类把 ``self.config`` 记为 ``TreePLConfig``，适配器后续只读子配置。
        # ``isinstance`` 与原先的 ``hasattr(config, "treepl")`` 对这两种入参等价：
        # ``ToolConfig`` 必有 ``treepl``，``TreePLConfig`` 必无。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            method_config = config.treepl
            if common is None:
                common = config.common
        else:
            method_config = config
        super().__init__(
            method_config, output_dir, software_paths, common_config=common
        )

        self.logger = get_logger()
        self._input_files: Dict[str, Path] = {}
        self._calibrations: List[CalibrationPoint] = []
        # 名称映射与（可能已短名化的）叶名：都在 prepare_inputs 里写，但
        # parse_results / CV 判定在“未调过 prepare_inputs”时也会读到它们，
        # 因此这里就把“可以为空”这件事说在类型里。
        self._name_mapping: Optional[Dict[str, str]] = None
        self._tree_tip_names: List[str] = []
        self._optimal_smooth: Optional[float] = None
        self._prime_result: Optional[PrimeResult] = None
        self._cv_results: List[CVResult] = []
        self._tree_modified: bool = False  # 标记是否修改了分支长度
        self._original_tree: Optional[PhylogeneticTree] = None
        self._actual_multiplier: float = 1.0  # 实际使用的分支长度乘数
        self._warnings: List[Any] = (
            []
        )  # 收集需要写入 DatingResult/runtime_metadata 的警告
        # 约束落盘台账（B-12 产物层口径）：written/skipped 都按校准名记
        self._constraint_ledger: Dict[str, List[str]] = {}

        # 设置随机种子
        if self.common_config.seed is None:
            self.common_config.seed = secrets.randbelow(1_000_000) + 1

    @property
    def method_name(self) -> str:
        return "treepl"

    def _get_executable(self) -> str:
        """获取 treePL 可执行文件路径，优先使用 software_paths 自定义路径。"""
        path = self.get_software_path("treepl_bin")
        return path if path else self.config.treepl_bin

    def _require_original_tree(self) -> PhylogeneticTree:
        """取 prepare_inputs 存下的输入树；未准备过就报错。

        这一步是调用顺序的前置条件（先 prepare_inputs 再 execute），不是数据
        质量问题：拿不到输入树时必须把它当成编程错误抬上去，而不是退回一个
        假的叶节点数（那会把"小树跳过 CV"的判定错算）。
        """
        if self._original_tree is None:
            raise ExecutionError(
                "treePL: the input tree is not available; prepare_inputs() must run "
                "before the CV/prime decision reads it"
            )
        return self._original_tree

    def validate_environment(self) -> bool:
        """检测 treePL 是否可用"""
        from ..infrastructure.software_dependencies import SoftwareDependencyManager

        runner = ProcessRunner()
        executable = self._get_executable()
        available = runner.check_executable(executable)

        if available:
            version = runner.get_version(executable, "--version")
            # treePL 不支持 --version；尝试从二进制中读取内嵌版本字符串
            if not version:
                version = self._extract_treepl_version(executable)
            self.software_version = version or "available (no --version support)"
            self.logger.info(f"treePL detected: {version or 'available'}")
        else:
            self.logger.warning(f"treePL not found: {executable}")
            install_cmd = SoftwareDependencyManager.get_install_command("treepl")
            self.logger.warning(f"Installation suggestions:\n{install_cmd}")

        return available

    def _extract_treepl_version(self, executable: str) -> Optional[str]:
        """从 treePL 二进制文件中提取版本字符串。

        treePL 通常不响应 --version，但二进制中常嵌入类似
        ``treePL version 1.0`` 的字符串。
        """
        import shutil

        path = shutil.which(executable)
        if not path:
            return None
        try:
            with open(path, "rb") as f:
                # 只读取前 1 MB，避免大文件耗时
                data = f.read(1_048_576)
            # 在可打印 ASCII 范围内搜索版本模式
            text = data.decode("ascii", errors="ignore")
            match = re.search(r"treePL\s+version\s+[\d.]+", text, re.IGNORECASE)
            if match:
                return match.group(0)
        except OSError as e:
            # 只影响版本溯源，不影响结果；但按 C-8 的要求不得静默吞掉
            self.logger.debug(f"Could not read treePL binary {path}: {e}")
        return None

    def _resolve_version_for_provenance(self, runner: ProcessRunner) -> str:
        """取 treePL 版本串用于执行溯源，失败时返回 "unknown" 并留痕（C-8）。

        过去 prime / CV / 最终分析三处各自写着
        `try: runner.get_version(...) except Exception: pass`，
        同一个探测失败在三条路径上都无声无息。
        """
        try:
            version_result = runner.get_version(self._get_executable(), "--version")
            if version_result:
                return version_result
        except Exception as e:
            self.logger.debug(f"treePL version probe failed: {e}")
        return "unknown"

    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
        生成 treePL 配置文件

        生成三个阶段的配置文件：
        1. Prime配置文件
        2. CV配置文件
        3. 最终定年配置文件
        """
        self._calibrations = calibrations
        self._original_tree = tree
        self._constraint_ledger = {}

        # treePL 对 numsites=0 的处理存在数值缺陷：当 numsites 为 0 时，
        # process_initial_branch_lengths 中的 min/numsites 会变成 inf，
        # 导致所有分支被设为 inf，进而 start_rate/initial calc 为 nan。
        # 如果用户未手动设置 numsites，自动从序列比对中检测并填入位点数。
        if (
            self.config.numsites == 0
            and alignment_path is not None
            and alignment_path.exists()
        ):
            extractor = AlignmentMetadataExtractor()
            metadata = extractor.extract(alignment_path)
            self.config.numsites = metadata.sequence_length
            self.logger.info(
                f"treePL: auto-setting numsites from alignment: {self.config.numsites}"
            )
        # 最后的兜底：确保 numsites > 0，避免 treePL 内部除零
        if self.config.numsites == 0:
            self.config.numsites = 1

        # 检查并处理 tiny branch
        processed_tree = self._check_and_fix_tiny_branches(tree)

        # 建立统一名称映射（树与 mrca/约束文件用同一套短名）。
        # treePL 对 tip 名长度/字符敏感，过长的 GTDB 全名会触发
        # "Unexisting tree file or Malformed newick tree structure" 等问题。
        name_mapping = self._build_name_mapping(processed_tree)
        if name_mapping is not None:
            processed_tree = processed_tree.with_renamed_leaves(name_mapping)
            self._name_mapping = name_mapping
            # 注意：不修改 calibrations 里的 mrca_leaf_pair（保持原始全名）。
            # _sanitize_name 已经会使用 _name_mapping 把全名映射为短名，
            # 用于生成 treePL 配置文件；parse_results 阶段则需要全名来定位
            # MRCA 并生成标准化输出。
        else:
            self._name_mapping = None

        # 保存输入树文件
        tree_file = self.work_dir / "input_tree.tre"
        processed_tree.without_internal_labels().write(tree_file)

        # 记录（可能被短名映射后的）全部叶节点名，供根节点校准生成 mrca 定义。
        # 树根 = 全部叶节点的 MRCA，故可用全部叶节点定义根节点的 mrca。
        self._tree_tip_names = list(processed_tree.tip_names)

        # 生成Prime配置文件
        prime_config_path = self._generate_prime_config(tree_file)

        # 生成CV配置文件
        cv_config_path = self._generate_cv_config(tree_file)

        self._input_files = {
            "tree": tree_file,
            "prime_config": prime_config_path,
            "cv_config": cv_config_path,
        }

        self.logger.info(f"Generated treePL configs in: {self.work_dir}")

        return self._input_files

    def _check_and_fix_tiny_branches(self, tree: PhylogeneticTree) -> PhylogeneticTree:
        """
        检查并处理 tiny branch

        根据教程 Step 4.5:
        如果最短分支长度小于 treePL 容忍的最小值，需要将所有分支长度
        乘以一个因子，使最短分支长度大于阈值。
        """
        if not self.config.check_tiny_branches:
            return tree

        # 获取所有分支长度
        branch_lengths = self._extract_branch_lengths(tree)

        if not branch_lengths:
            return tree

        min_bl = min(branch_lengths)
        self.logger.info(
            f"Tree branch length stats: min={min_bl:.2e}, max={max(branch_lengths):.2e}"
        )

        # 检查是否需要调整
        if min_bl < self.config.min_branch_length_threshold:
            # 计算需要的乘数
            multiplier = self.config.branch_length_multiplier
            required_multiplier = (
                self.config.min_branch_length_threshold * 10
            ) / min_bl

            if required_multiplier > multiplier:
                multiplier = required_multiplier

            # 保存实际使用的乘数（用于后续回退调整）
            self._actual_multiplier = multiplier

            self.logger.warning(
                f"Detected tiny branches (min={min_bl:.2e}). "
                f"Multiplying all branch lengths by {multiplier:.0f}"
            )

            # 修改分支长度
            modified_tree = self._multiply_branch_lengths(tree, multiplier)

            # 检查是否成功修改（ETE3 可用时才会真正修改）
            if modified_tree is not tree:
                self._tree_modified = True
                self._actual_multiplier = multiplier
            else:
                # ETE3 不可用，无法修改分支长度
                self._tree_modified = False
                self._actual_multiplier = 1.0
                self.logger.warning(
                    "ETE3 not available, cannot modify branch lengths. "
                    "Tiny branch handling skipped."
                )

            # 保存原始树
            original_tree_file = self.work_dir / "original_tree.tre"
            tree.write(original_tree_file)
            self.logger.info(f"Original tree saved to: {original_tree_file}")

            # 验证修改后的树
            new_lengths = self._extract_branch_lengths(modified_tree)
            self.logger.info(
                f"After modification: min={min(new_lengths):.2e}, "
                f"max={max(new_lengths):.2e}"
            )

            return modified_tree

        return tree

    def _extract_branch_lengths(self, tree: PhylogeneticTree) -> List[float]:
        """提取树中所有分支长度"""
        lengths = []
        try:
            from ete3 import Tree

            ete_tree = Tree(strip_leading_newick_comments(tree.newick), format=1)

            for node in ete_tree.traverse():
                if node.dist is not None and node.dist > 0:
                    lengths.append(node.dist)

            return lengths
        except ImportError:
            self.logger.warning("ETE3 not available for branch length extraction")
            return []

    def _multiply_branch_lengths(
        self, tree: PhylogeneticTree, multiplier: float
    ) -> PhylogeneticTree:
        """将所有分支长度乘以指定因子"""
        try:
            from ete3 import Tree

            # 复制树
            ete_tree = Tree(strip_leading_newick_comments(tree.newick), format=1)

            # 修改所有分支长度
            for node in ete_tree.traverse():
                if node.dist is not None:
                    node.dist *= multiplier

            # 写出新的 Newick
            new_newick = ete_tree.write(format=5)  # format=5: 保留分支长度

            return PhylogeneticTree.from_newick(new_newick)

        except ImportError:
            self.logger.warning("ETE3 not available, cannot modify branch lengths")
            return tree

    def _write_common_params(self, f: TextIO, tree_file: Path) -> None:
        """写入通用参数（所有配置文件共享）"""
        # 树文件
        f.write(f"treefile = {tree_file}\n")

        # 树处理参数
        if self.config.numsites > 0:
            f.write(f"numsites = {self.config.numsites}\n")
        if self.config.scale != 1.0:
            f.write(f"scale = {self.config.scale}\n")
        if self.config.collapse:
            f.write("collapse\n")
        if self.config.set1:
            f.write("set1\n")

        # 线程数（使用 common 配置）
        f.write(f"nthreads = {self.common_config.nthreads}\n")

        # 随机种子（使用 common 配置）
        if self.common_config.seed is not None:
            f.write(f"seed = {self.common_config.seed}\n")

        # 优化参数
        if self.config.opt is not None:
            f.write(f"opt = {self.config.opt}\n")
        if self.config.optad is not None:
            f.write(f"optad = {self.config.optad}\n")
        if self.config.optcvad is not None:
            f.write(f"optcvad = {self.config.optcvad}\n")

        # 详细优化参数
        if self.config.moredetail:
            f.write("moredetail\n")
        if self.config.moredetailad:
            f.write("moredetailad\n")
        if self.config.moredetailcvad:
            f.write("moredetailcvad\n")

        # 优化器参数
        if self.config.calcgrad:
            f.write("calcgrad\n")
        if self.config.thorough:
            f.write("thorough\n")
        if self.common_config.verbose:
            f.write("verbose\n")
        if self.config.paramverbose:
            f.write("paramverbose\n")

        # 收敛容差
        if self.config.ftol is not None:
            f.write(f"ftol = {self.config.ftol}\n")
        if self.config.xtol is not None:
            f.write(f"xtol = {self.config.xtol}\n")

        # 模拟退火参数
        if self.config.lftemp is not None:
            f.write(f"lftemp = {self.config.lftemp}\n")
        if self.config.pltemp is not None:
            f.write(f"pltemp = {self.config.pltemp}\n")
        if self.config.lfcool is not None:
            f.write(f"lfcool = {self.config.lfcool}\n")
        if self.config.plcool is not None:
            f.write(f"plcool = {self.config.plcool}\n")
        if self.config.lfstoptemp is not None:
            f.write(f"lfstoptemp = {self.config.lfstoptemp}\n")
        if self.config.plstoptemp is not None:
            f.write(f"plstoptemp = {self.config.plstoptemp}\n")
        if self.config.lfrtstep is not None:
            f.write(f"lfrtstep = {self.config.lfrtstep}\n")
        if self.config.lfdtstep is not None:
            f.write(f"lfdtstep = {self.config.lfdtstep}\n")
        if self.config.plrtstep is not None:
            f.write(f"plrtstep = {self.config.plrtstep}\n")
        if self.config.pldtstep is not None:
            f.write(f"pldtstep = {self.config.pldtstep}\n")

        # 迭代次数参数
        if self.config.lfiter is not None:
            f.write(f"lfiter = {self.config.lfiter}\n")
        if self.config.pliter is not None:
            f.write(f"pliter = {self.config.pliter}\n")
        if self.config.lfsimaniter is not None:
            f.write(f"lfsimaniter = {self.config.lfsimaniter}\n")
        if self.config.plsimaniter is not None:
            f.write(f"plsimaniter = {self.config.plsimaniter}\n")

        # 其他参数
        if self.config.checkconstraints:
            f.write("checkconstraints\n")
        if self.config.mapspace:
            f.write("mapspace\n")
        if self.config.log_pen:
            f.write("log_pen\n")
        if self.config.sample != 1.0:
            f.write(f"sample = {self.config.sample}\n")

        # 输入文件参数
        if self.config.ind8s is not None:
            f.write(f"ind8s = {self.config.ind8s}\n")
        if self.config.inr8s is not None:
            f.write(f"inr8s = {self.config.inr8s}\n")

    def _write_mrca_and_constraints(self, f: TextIO) -> None:
        """写入MRCA定义和约束

        优先使用 ``mrca_leaf_pair``（两个叶节点名），若为空则回退到
        ``resolved_taxa``。这样即使校准加载器未填充 ``resolved_taxa``，
        只要命令行 YAML 里写了 ``mrca_pair`` 就能正确生成 MRCA。

        对于 ``is_root_node`` 的根节点校准（无 mrca_leaf_pair），用树的全部
        叶节点定义 mrca（根 = 全部叶节点的最近共同祖先），否则 treePL 会因
        只写 ``max = Root X`` 而无对应 mrca 定义而报错。

        每条没能落盘的校准都会被记名并警告一次（本方法在 prime/CV/final
        三阶段被反复调用，用 `_warnings_set` 去重），并由
        ``_verify_written_config`` 在写盘后按产物复核。
        """
        self._ensure_constraint_ledger()
        for cal in self._calibrations:
            # 优先 mrca_leaf_pair，回退 resolved_taxa
            taxa = (
                list(cal.mrca_leaf_pair)
                if cal.mrca_leaf_pair
                else list(cal.resolved_taxa or [])
            )

            # 根节点校准：无叶节点对，用全部叶节点定义 mrca
            if not taxa and cal.is_root_node:
                tips = getattr(self, "_tree_tip_names", None)
                if tips:
                    taxa = list(tips)
                    self.logger.info(
                        f"treePL: defining root calibration '{cal.name}' as mrca of all {len(tips)} tips"
                    )

            if not taxa:
                self._note_unwritten_calibration(cal)
                continue

            if cal.age_constraint is None:
                # 与"定位不到叶对"同一条口径：没有年龄界的校准点写不出
                # minmax/fixed 行，记账后跳过，而不是在 None 上调
                # to_software_format 碰 AttributeError。
                self._note_unwritten_calibration(cal)
                continue

            taxa_str = " ".join(self._sanitize_name(t) for t in taxa)
            f.write(f"mrca = {cal.name} {taxa_str}\n")

            self._write_treepl_constraint(f, cal)
            self._constraint_ledger["written"].append(cal.name)

    def _ensure_constraint_ledger(self) -> None:
        if "written" not in self._constraint_ledger:
            self._constraint_ledger["written"] = []
        if "skipped" not in self._constraint_ledger:
            self._constraint_ledger["skipped"] = []

    def _record_degradation_warning(self, msg: str) -> None:
        """记录一条需要随结果交出去的警告，并去重。

        本适配器的写盘方法会在 prime / CV / final 三个阶段被反复调用，
        不去重就会让同一个问题在结果里出现三遍。
        """
        if not hasattr(self, "_warnings_set"):
            self._warnings_set: Set[str] = set()
        if msg in self._warnings_set:
            return
        self._warnings_set.add(msg)
        self._warnings.append(SemanticDegradationWarning(msg))
        self.logger.warning(msg, category=SemanticDegradationWarning)

    def _note_unwritten_calibration(self, cal: CalibrationPoint) -> None:
        """一条校准既无 mrca_leaf_pair 也无 resolved_taxa：不能悄悄消失。"""
        self._ensure_constraint_ledger()
        if cal.name in self._constraint_ledger["skipped"]:
            return
        self._constraint_ledger["skipped"].append(cal.name)
        msg = (
            f"Calibration '{cal.name}' was NOT written to the treePL configuration "
            "file: it has neither mrca_leaf_pair nor resolved_taxa, so no mrca can "
            "be defined for it. It has no effect on the dating result."
        )
        self._record_degradation_warning(msg)

    def _verify_written_config(self, config_path: Path) -> None:
        """回读真正交给 treePL 的配置文件，确认约束确实落盘（B-12 的产物层口径）。

        treePL 没有 PATHd8 那种"根=1"的失败模式，但它同样会在配置里一条
        mrca/约束都没有时安静地拟合一棵不受任何校准约束的树；而且 treePL 在
        多种错误路径上打印信息后 `exit(0)`（审阅报告 §六 末"共性风险"），
        所以判据必须是**写出去的产物**，而不是内存里的校准列表。
        """
        try:
            with open(config_path, "r") as f:
                text = f.read()
        except OSError as e:
            raise ExecutionError(
                f"treePL configuration file {config_path} could not be read back "
                f"for verification: {e}"
            )

        mrca_lines = len(
            re.findall(r"^\s*mrca\s*=", text, re.IGNORECASE | re.MULTILINE)
        )
        bound_lines = len(
            re.findall(
                r"^\s*(?:min|max)\s*=\s*\S+\s+[-+0-9.eE]+",
                text,
                re.IGNORECASE | re.MULTILINE,
            )
        )

        if mrca_lines == 0 or bound_lines == 0:
            skipped = self._constraint_ledger.get("skipped") or []
            raise CalibrationError(
                f"treePL configuration '{config_path.name}' carries no usable "
                f"calibration (mrca definitions: {mrca_lines}, min/max bound lines: "
                f"{bound_lines}) even though {len(self._calibrations)} calibration(s) "
                "were supplied. treePL would fit an unconstrained tree here and can "
                "still exit 0, so this is refused up front. Check that every "
                "calibration resolves to an mrca present in the tree "
                "(mrca_leaf_pair / resolved_taxa)."
                + (
                    f" Not written: {', '.join(sorted(set(skipped)))}"
                    if skipped
                    else ""
                )
            )

        if len(self._calibrations) > len(
            set(self._constraint_ledger.get("written") or [])
        ):
            msg = (
                f"treePL input reconciliation: {len(self._calibrations)} "
                f"calibration(s) supplied, "
                f"{len(set(self._constraint_ledger.get('written') or []))} of them "
                "reached the configuration file."
            )
            self._record_degradation_warning(msg)

    def _generate_prime_config(self, tree_file: Path) -> Path:
        """生成Prime步骤配置文件"""
        prime_config_path = self.work_dir / "treepl_prime.conf"

        with safe_writer(prime_config_path) as f:
            # 写入通用参数
            self._write_common_params(f, tree_file)

            # Prime命令 - 获取优化参数
            f.write("prime\n")

            # MRCA定义和约束（Prime也需要）
            self._write_mrca_and_constraints(f)

        self._verify_written_config(prime_config_path)
        self.logger.debug(f"Generated prime config: {prime_config_path}")
        return prime_config_path

    def _generate_cv_config(self, tree_file: Path) -> Path:
        """生成CV步骤配置文件"""
        cv_config_path = self.work_dir / "treepl_cv.conf"

        with safe_writer(cv_config_path) as f:
            # 写入通用参数
            self._write_common_params(f, tree_file)

            # CV输出文件
            f.write(f"cvoutfile = {self.work_dir / self.config.cvoutfile}\n")

            # CV类型：randomcv（推荐）或 cv
            if self.config.randomcv:
                f.write("randomcv\n")
            elif self.config.cv:
                f.write("cv\n")

            # CV参数
            f.write(f"cvstart = {self.config.cvstart}\n")
            f.write(f"cvstop = {self.config.cvstop}\n")
            f.write(f"cvmultstep = {self.config.cvmultstep}\n")

            # randomcv特有参数
            if self.config.randomcv:
                f.write(f"cviter = {self.config.cviter}\n")
                if self.config.cvsimaniter > 0:
                    f.write(f"cvsimaniter = {self.config.cvsimaniter}\n")
                if self.config.randomcvsamp != 0.1:
                    f.write(f"randomcvsamp = {self.config.randomcvsamp}\n")

            # MRCA定义和约束
            self._write_mrca_and_constraints(f)

        self._verify_written_config(cv_config_path)
        self.logger.debug(f"Generated CV config: {cv_config_path}")
        return cv_config_path

    def _generate_final_config(self, tree_file: Path, smooth: float) -> Path:
        """生成最终定年配置文件"""
        final_config_path = self.work_dir / "treepl_final.conf"

        with safe_writer(final_config_path) as f:
            # 写入通用参数
            self._write_common_params(f, tree_file)

            # 输出文件
            f.write(f"outfile = {self.work_dir / self.config.outfile}\n")

            # 平滑参数
            f.write(f"smooth = {smooth}\n")

            # MRCA定义和约束
            self._write_mrca_and_constraints(f)

        self._verify_written_config(final_config_path)
        self.logger.debug(f"Generated final config: {final_config_path}")
        return final_config_path

    def _build_name_mapping(self, tree: PhylogeneticTree) -> Optional[Dict[str, str]]:
        """建立全名 → 安全短名的映射，使树与 mrca/约束文件名称一致。

        treePL 对 tip 名长度和字符敏感，GTDB 全名（含 ``.`` ``-`` 和很长）会
        触发 "Malformed newick tree structure" 等问题。这里用 NameMappingManager
        的 ACCESSION 模式提取 ``GB_GCA_<accession>`` 部分（既短又唯一）。
        如果所有名本身已 ≤30 字符且 PHYLIP 安全，则无需映射。
        """
        from ..services.name_mapping import NameMappingManager, NameShortenMode

        mapping = NameMappingManager()
        any_changed = False
        for original in tree.tip_names:
            short = mapping.add_entry(original, method=NameShortenMode.ACCESSION)
            if short != original:
                any_changed = True
        return mapping.get_mapping_dict() if any_changed else None

    def _sanitize_name(self, name: str) -> str:
        """安全化名称，使其既 ≤30 字符又与树 tip 名一致。

        treePL 内部会把 tip 名截断到 30 字符，所以这里的截断必须和树里的
        名同形。策略：用映射表把全名替换成 ACCESSION 短名（保证 ≤18 字符）；
        若没有映射则回退到传统的正则截断。
        """
        mapping = self._name_mapping
        if mapping and name in mapping:
            return mapping[name]
        safe = re.sub(r"[^A-Za-z0-9_]", "_", name)
        return safe[:30]

    def _write_treepl_constraint(self, f: TextIO, cal: CalibrationPoint) -> None:
        """将校准约束写入 treePL 配置文件

        统一使用约束对象的 to_software_format('treepl') 生成，仅在需要时补充
        适配器级降级警告。
        """
        constraint = cal.age_constraint
        if constraint is None:
            # 调用方（_write_mrca_and_constraints）已在写 mrca 行前拦掉；
            # 走到这里说明调用顺序被破坏，报错比拿 None 去渲染更可诊断。
            raise CalibrationError(
                f"treePL: calibration '{cal.name}' has no age_constraint; nothing can be "
                "written for it"
            )

        # 直接传入校准点名，模型内部完成占位符渲染。
        formatted = constraint.to_software_format("treepl", name=cal.name)
        f.write(f"{formatted}\n")

        # 为需要降级的约束类型补充适配器级警告（to_software_format 已发出模型级警告）。
        # 由于本方法会被 prime/cv/final 三阶段反复调用，_record_degradation_warning 去重。
        def _add_warning(msg: str) -> None:
            self._record_degradation_warning(msg)

        if isinstance(constraint, SoftLowerBoundConstraint):
            _add_warning(
                f"Soft lower bound degraded to hard min for treePL: {cal.name}"
            )
        elif isinstance(constraint, MaximumAgeConstraint):
            _add_warning(f"Maximum age degraded to hard max for treePL: {cal.name}")
        elif isinstance(
            constraint,
            (
                SoftBoundsConstraint,
                GammaPriorConstraint,
                SkewNormalConstraint,
                SkewTConstraint,
            ),
        ):
            _add_warning(
                f"{type(constraint).__name__} degraded to uniform range for treePL: {cal.name}"
            )

    def execute(self) -> bool:
        """
        执行 treePL 分析（三阶段流程）

        阶段 1: Prime - 获取最佳优化参数(opt/optad)
        阶段 2: CV分析 - 获取最佳平滑参数
        阶段 3: 最终定年 - 使用优化参数和平滑参数
        """
        tree_file = self._input_files["tree"]

        # 阶段 1: Prime（如果启用且参数未手动设置）
        if self.config.run_prime and (
            self.config.opt is None or self.config.optad is None
        ):
            self.logger.step(1, "Prime - Finding optimal optimization parameters")
            prime_result = self._run_prime_with_replicates()

            if prime_result.success:
                self._prime_result = prime_result
                self.config.opt = prime_result.opt
                self.config.optad = prime_result.optad
                self.logger.info(
                    f"Prime result: opt={prime_result.opt}, optad={prime_result.optad}"
                )

                # 重新生成CV配置（包含优化参数）
                self._generate_cv_config(tree_file)
            else:
                # Prime 全部 replicate 失败通常是配置/输入缺陷（确定性错误）。
                # 根据 treePL 输出的具体信息区分根因，给出更准确的建议。
                msg_lower = prime_result.message.lower()
                stdout_lower = (prime_result.stdout or "").lower()
                combined = msg_lower + " " + stdout_lower

                if (
                    "feasible start rates" in combined
                    or "initial calc: nan" in combined
                ):
                    # 优化器无法找到可行初始速率（常见于叶节点少或约束稀疏的小树）
                    suggestion = (
                        "This is a treePL optimizer initialization failure: it could not "
                        "find feasible starting rates/dates for the given tree and calibrations. "
                        "Try adding more calibration points, using a less sparse constraint set, "
                        "or choose another dating method (e.g. r8s/lsd2)."
                    )
                elif "could not open file" in combined or "incomplete" in combined:
                    # 历史问题：mrca 定义缺失导致配置文件中存在空文件名
                    suggestion = (
                        "This typically indicates an incomplete calibration mrca definition. "
                        "Please check that every calibration has a valid mrca_pair and that the "
                        "tip names match the tree labels."
                    )
                else:
                    suggestion = (
                        "Please check your calibrations and tree, add more calibration "
                        "points, or use a different dating method."
                    )

                raise ExecutionError(
                    f"treePL Prime step failed: {prime_result.message}. "
                    "All prime replicates failed, which indicates a deterministic "
                    "configuration or input problem. CV and final analysis are skipped "
                    "to avoid redundant failures. "
                    f"{suggestion}"
                )
        else:
            self.logger.info("Skipping prime step (disabled or manually configured)")
            if self.config.opt is None:
                self.config.opt = 1
            if self.config.optad is None:
                self.config.optad = 1

        # 阶段 2: 交叉验证（可重复运行以检查稳定性）
        cv_type_str = "randomcv" if self.config.randomcv else "cv"
        self.logger.step(2, f"Cross-validation ({cv_type_str}) for smoothing parameter")

        # 当所有校准都是固定年龄时，randomcv/cv 无法 hold-out，直接采用 initial_smooth。
        all_fixed = all(
            isinstance(cal.age_constraint, FixedAgeConstraint)
            for cal in self._calibrations
        )
        # 对于叶节点很少的小树，交叉验证同样缺乏足够的 hold-out 折叠，
        # 直接采用 initial_smooth 可避免无意义的 CV 失败警告。
        n_tips = len(self._tree_tip_names or self._require_original_tree().tip_names)
        small_tree = n_tips <= 10
        if all_fixed or small_tree:
            reason = (
                "all calibrations are fixed ages"
                if all_fixed
                else f"small tree (n_tips={n_tips})"
            )
            self.logger.info(
                f"{reason}; skipping CV and using initial_smooth = {self.config.initial_smooth}"
            )
            self._optimal_smooth = self.config.initial_smooth
            self.logger.step(3, "Final dating analysis with optimal parameters")
            self._run_final_analysis(self._optimal_smooth)
            return True

        self.logger.info(
            f"Running {self.config.cv_replicates} CV replicates for stability check"
        )

        # 如果需要，自动调整 opt/optad
        if self.config.auto_adjust_opt:
            self._auto_adjust_opt_params(tree_file)

        cv_smooth_values = []
        for i in range(self.config.cv_replicates):
            if self.config.cv_replicates > 1:
                self.logger.info(f"CV replicate {i+1}/{self.config.cv_replicates}")

            optimal_smooth = self._run_cross_validation()

            if optimal_smooth is not None:
                cv_smooth_values.append(optimal_smooth)
                self._cv_results.append(
                    CVResult(
                        smoothing=optimal_smooth,
                        chisq=0.0,  # 稍后从文件解析
                        success=True,
                    )
                )

        if cv_smooth_values:
            # 根据教程建议：如果预期有大量速率异质性，使用最低值
            self._optimal_smooth = min(cv_smooth_values)
            self.logger.info(f"CV results: {cv_smooth_values}")
            self.logger.info(
                f"Selected optimal smoothing (lowest for rate heterogeneity): {self._optimal_smooth}"
            )

            if len(cv_smooth_values) > 1:
                import statistics

                mean_smooth = statistics.mean(cv_smooth_values)
                std_smooth = (
                    statistics.stdev(cv_smooth_values)
                    if len(cv_smooth_values) > 1
                    else 0
                )
                self.logger.info(
                    f"CV stability: mean={mean_smooth:.6f}, std={std_smooth:.6f}"
                )

            # 自适应调整CV参数
            if self.config.adaptive_cv:
                self._adaptive_cv_adjustment(cv_smooth_values)
        else:
            self.logger.warning(
                f"All CV runs failed, using initial smooth: {self.config.initial_smooth}",
                category=ConvergenceWarning,
            )
            self._optimal_smooth = self.config.initial_smooth

        # 阶段 3: 最终分析
        self.logger.step(3, "Final dating analysis with optimal parameters")
        self._run_final_analysis(self._optimal_smooth)

        return True

    def _run_prime_with_replicates(self) -> PrimeResult:
        """
        多次运行Prime步骤，选择最低的opt和optad

        根据教程: "Repeat the priming analysis a few time, and select the lines
        that have the lowest values for opt and optad"
        """
        all_results: List[PrimeResult] = []

        for i in range(self.config.prime_replicates):
            self.logger.info(f"Prime replicate {i+1}/{self.config.prime_replicates}")
            result = self._run_prime()

            if result.success:
                all_results.append(result)
                self.logger.info(
                    f"  Prime {i+1}: opt={result.opt}, optad={result.optad}"
                )

        if not all_results:
            # 保留最后一次失败的 stdout，帮助上层区分是 mrca 配置问题还是
            # treePL 数值初始化问题。
            last_stdout = ""
            if hasattr(self, "_last_prime_stdout"):
                last_stdout = self._last_prime_stdout
            return PrimeResult(
                opt=1,
                optad=1,
                success=False,
                message="All prime replicates failed",
                stdout=last_stdout,
            )

        # 选择最低的opt和optad
        min_opt = min(r.opt for r in all_results)
        min_optad = min(r.optad for r in all_results)

        self.logger.info(
            f"Selected lowest values from {len(all_results)} replicates: opt={min_opt}, optad={min_optad}"
        )

        return PrimeResult(opt=min_opt, optad=min_optad, success=True)

    def _run_prime(self) -> PrimeResult:
        """
        运行Prime步骤获取优化参数

        解析输出中的PLACE THE LINES BELOW IN THE CONFIGURATION FILE部分
        提取opt和optad值
        """
        prime_config = self._input_files["prime_config"]

        cmd = [self._get_executable(), str(prime_config)]

        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 7200
        )

        version_info = self._resolve_version_for_provenance(runner)

        self.logger.execution_provenance(
            tool_name="treePL",
            command=" ".join(cmd),
            version=version_info,
            environment={
                "treePL_version": version_info,
                "opt": str(self.config.opt),
                "optad": str(self.config.optad) if self.config.optad else "auto",
                "work_dir": str(self.work_dir),
            },
            artifacts={
                "input_tree": str(self._input_files["tree"]),
                "prime_config": str(prime_config),
            },
        )

        result = runner.run(cmd)

        if result.returncode != 0:
            self._last_prime_stdout = result.stdout
            return PrimeResult(
                opt=1,
                optad=1,
                success=False,
                message=f"Prime failed with returncode {result.returncode}: {result.stderr}",
                stdout=result.stdout,
            )

        # 解析Prime输出，寻找优化参数
        parsed = self._parse_prime_output(result.stdout)
        if not parsed.success:
            self._last_prime_stdout = result.stdout
        return parsed

    def _parse_prime_output(self, output: str) -> PrimeResult:
        """
        解析Prime输出，提取opt和optad参数

        寻找格式：
        PLACE THE LINES BELOW IN THE CONFIGURATION FILE
        opt = X
        optad = Y
        """
        # 查找配置建议部分
        config_section = re.search(
            r"PLACE THE LINES BELOW IN THE CONFIGURATION FILE\s*\n(.*?)\n\n",
            output,
            re.IGNORECASE | re.DOTALL,
        )

        if config_section:
            section = config_section.group(1)

            # 提取opt
            opt_match = re.search(r"opt\s*=\s*(\d+)", section, re.IGNORECASE)
            optad_match = re.search(r"optad\s*=\s*(\d+)", section, re.IGNORECASE)

            if opt_match and optad_match:
                opt = int(opt_match.group(1))
                optad = int(optad_match.group(1))
                return PrimeResult(opt=opt, optad=optad, success=True)

        # 备选：直接在完整输出中搜索
        opt_match = re.search(r"opt\s*=\s*(\d+)", output, re.IGNORECASE)
        optad_match = re.search(r"optad\s*=\s*(\d+)", output, re.IGNORECASE)

        if opt_match and optad_match:
            opt = int(opt_match.group(1))
            optad = int(optad_match.group(1))
            return PrimeResult(opt=opt, optad=optad, success=True)

        return PrimeResult(
            opt=1,
            optad=1,
            success=False,
            message="Could not parse opt/optad from prime output",
            stdout=output,
        )

    def _auto_adjust_opt_params(self, tree_file: Path) -> None:
        """
        自动调整opt/optad以消除警告消息

        根据教程 Step 5.5:
        当出现 "(might want to try a different opt=VALUE)" 或
        "(might want to try a different optad=VALUE)" 消息时，
        需要增加相应的值直到消息消失。
        """
        max_attempts = self.config.max_opt_adjust
        attempt = 0

        while attempt < max_attempts:
            # 运行CV检查是否有警告
            cv_config = self._generate_cv_config(tree_file)
            cmd = [self._get_executable(), str(cv_config)]

            runner = ProcessRunner(
                cwd=self.work_dir, timeout=self.common_config.timeout or 3600
            )
            result = runner.run(cmd)

            if result.returncode != 0:
                self.logger.error(f"CV check failed: {result.stderr}")
                break

            # 检查警告消息
            output = result.stdout + result.stderr
            opt_warning = "might want to try a different opt=" in output.lower()
            optad_warning = "might want to try a different optad=" in output.lower()

            if not opt_warning and not optad_warning:
                self.logger.info("No opt/optad warnings detected")
                break

            # 增加相应的值。opt/optad 只有在"当前值已知"时才能递：prime 阶段
            # （或用户显式配置）会给出数值；仍为 None 说明用的一直是 treePL 自己
            # 的内置默认，我们无从知道该从几开始，只能把它说清楚并停止这一项
            # 调整（旧写法 ``None + 1`` 会直接 TypeError）。
            if opt_warning:
                if self.config.opt is None:
                    self.logger.warning(
                        "treePL: the output asks for a different 'opt', but opt is "
                        "unset (no prime result and no user value), so it cannot be "
                        "increased automatically. Set software.treepl.opt to hand "
                        "the automatic adjustment a starting value."
                    )
                else:
                    self.config.opt += 1
                    self.logger.info(
                        f"Detected opt warning, increasing opt to {self.config.opt}"
                    )

            if optad_warning:
                if self.config.optad is None:
                    self.logger.warning(
                        "treePL: the output asks for a different 'optad', but optad "
                        "is unset (no prime result and no user value), so it cannot "
                        "be increased automatically. Set software.treepl.optad to "
                        "hand the automatic adjustment a starting value."
                    )
                else:
                    self.config.optad += 1
                    self.logger.info(
                        f"Detected optad warning, increasing optad to "
                        f"{self.config.optad}"
                    )

            attempt += 1

            # 重新生成配置
            self._generate_cv_config(tree_file)

        if attempt >= max_attempts:
            self.logger.warning(
                f"Reached max adjustment attempts ({max_attempts}). "
                f"Using opt={self.config.opt}, optad={self.config.optad}"
            )
        else:
            self.logger.info(
                f"Opt/optad adjustment complete: opt={self.config.opt}, optad={self.config.optad}"
            )

    def _adaptive_cv_adjustment(self, cv_smooth_values: List[float]) -> None:
        """
        自适应调整CV参数

        根据教程:
        如果cvstop值具有最低的chisq值，说明最佳平滑值可能低于当前cvstop，
        需要降低cvstop继续搜索。
        """
        if not self.config.adaptive_cv:
            return

        # 检查cvoutfile以获取chisq值
        cvoutfile_path = self.work_dir / self.config.cvoutfile
        if not cvoutfile_path.exists():
            return

        try:
            with open(cvoutfile_path, "r") as f:
                content = f.read()

            # 解析CV结果
            best_smooth, best_chisq = self._parse_cv_file_with_chisq(content)

            if best_smooth is None:
                return

            # 如果最佳平滑值接近cvstop，可能需要降低cvstop
            # 添加数值稳定性保护：使用 max(cvstop, cvstop_min) 避免除以极小数
            safe_cvstop = max(self.config.cvstop, self.config.cvstop_min, 1e-30)
            if abs(best_smooth - self.config.cvstop) / safe_cvstop < 0.1:
                new_cvstop = self.config.cvstop * 0.1

                if new_cvstop >= self.config.cvstop_min:
                    self.logger.info(
                        f"Best smoothing ({best_smooth}) close to cvstop ({self.config.cvstop}). "
                        f"Reducing cvstop to {new_cvstop} for better search"
                    )
                    self.config.cvstop = new_cvstop

                    # 重新运行一次CV
                    self.logger.info("Re-running CV with adjusted cvstop")
                    tree_file = self._input_files["tree"]
                    self._generate_cv_config(tree_file)

                    optimal_smooth = self._run_cross_validation()
                    if optimal_smooth is not None:
                        self._optimal_smooth = optimal_smooth
                        self.logger.info(
                            f"Updated optimal smoothing: {self._optimal_smooth}"
                        )

        except Exception as e:
            self.logger.warning(f"Adaptive CV adjustment failed: {e}")

    def _parse_cv_file_with_chisq(
        self, content: str
    ) -> Tuple[Optional[float], Optional[float]]:
        """解析CV文件，返回最佳平滑值和对应的chisq"""
        best_smooth = None
        best_chisq = float("inf")

        lines = content.strip().split("\n")
        data_started = False

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if "smoothing" in line.lower() and "chi" in line.lower():
                data_started = True
                continue

            if not data_started:
                continue

            parts = line.split()
            if len(parts) >= 2:
                try:
                    smooth = float(parts[0])
                    chisq = float(parts[1])

                    if chisq < best_chisq:
                        best_chisq = chisq
                        best_smooth = smooth
                except ValueError:
                    continue

        return best_smooth, best_chisq

    def _run_cross_validation(self) -> Optional[float]:
        """运行交叉验证"""
        cv_config = self._input_files["cv_config"]

        cmd = [self._get_executable(), str(cv_config)]

        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 7200
        )

        version_info = self._resolve_version_for_provenance(runner)

        self.logger.execution_provenance(
            tool_name="treePL",
            command=" ".join(cmd),
            version=version_info,
            environment={
                "treePL_version": version_info,
                "cv_type": "randomcv" if self.config.randomcv else "cv",
                "work_dir": str(self.work_dir),
            },
            artifacts={
                "input_tree": str(self._input_files["tree"]),
                "cv_config": str(cv_config),
            },
        )

        result = runner.run(cmd)

        if result.returncode != 0:
            self.logger.error(f"CV failed: {result.stderr}")
            return None

        # 首先尝试从cvoutfile解析
        cvoutfile_path = self.work_dir / self.config.cvoutfile
        if cvoutfile_path.exists():
            try:
                with open(cvoutfile_path, "r") as f:
                    smooth = self._parse_cv_file(f.read())
                if smooth is not None:
                    return smooth
            except Exception as e:
                self.logger.warning(f"Failed to parse cvoutfile: {e}")

        # 备选：从stdout解析
        smooth = self._parse_cv_output(result.stdout)
        return smooth

    def _parse_cv_file(self, content: str) -> Optional[float]:
        """
        解析CV输出文件，找到最优平滑参数

        格式示例：
        smoothing    chi-square
        1000.0       123.45
        100.0        89.12
        ...
        """
        best_smooth = None
        best_chisq = float("inf")

        lines = content.strip().split("\n")

        # 跳过标题行
        data_started = False
        for line in lines:
            line = line.strip()
            if not line:
                continue

            # 检测数据开始（跳过标题）
            if "smoothing" in line.lower() and "chi" in line.lower():
                data_started = True
                continue

            if not data_started:
                continue

            # 解析数据行
            parts = line.split()
            if len(parts) >= 2:
                try:
                    smooth = float(parts[0])
                    chisq = float(parts[1])

                    if chisq < best_chisq:
                        best_chisq = chisq
                        best_smooth = smooth
                except ValueError:
                    continue

        return best_smooth

    def _parse_cv_output(self, output: str) -> Optional[float]:
        """解析 CV 输出，找到最优平滑参数"""
        best_smooth = None
        best_chisq = float("inf")

        # treePL CV 输出格式示例：
        # smoothing = 100.0 chi-square = 123.45
        pattern = r"smoothing\s*=\s*([\d.eE+-]+).*?chi-square\s*=\s*([\d.eE+-]+)"

        for match in re.finditer(pattern, output, re.IGNORECASE):
            try:
                smooth = float(match.group(1))
                chisq = float(match.group(2))

                if chisq < best_chisq and chisq != float("inf"):
                    best_chisq = chisq
                    best_smooth = smooth
            except ValueError:
                continue

        return best_smooth

    def _run_final_analysis(self, smooth: float) -> None:
        """运行最终分析"""
        tree_file = self._input_files["tree"]

        final_config_path = self._generate_final_config(tree_file, smooth)

        cmd = [self._get_executable(), str(final_config_path)]

        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 7200
        )

        version_info = self._resolve_version_for_provenance(runner)

        self.logger.execution_provenance(
            tool_name="treePL",
            command=" ".join(cmd),
            version=version_info,
            environment={
                "treePL_version": version_info,
                "smoothing": str(smooth),
                "work_dir": str(self.work_dir),
            },
            artifacts={
                "input_tree": str(tree_file),
                "final_config": str(final_config_path),
                "output_dated_tree": str(self.work_dir / self.config.outfile),
            },
        )

        result = runner.run(cmd)

        if result.returncode != 0:
            raise ExecutionError(
                f"Final analysis failed (returncode={result.returncode}): {result.stderr}"
            )

        self._input_files["final_config"] = final_config_path
        self._input_files["dated_tree"] = self.work_dir / self.config.outfile

        self.logger.success("Final analysis completed")

    def parse_results(self) -> DatingResult:
        """解析 treePL 输出"""
        dated_tree_path = self._input_files.get("dated_tree")

        if not dated_tree_path or not dated_tree_path.exists():
            raise ResultParsingError("treePL dated tree not found")

        with open(dated_tree_path, "r") as f:
            dated_tree_newick = f.read().strip()

        # treePL 内部优化失败时会写出空文件并以 rc=0 退出。
        # 提前检查，给出更清晰的错误，避免 ETE3 解析空文件报
        # "Malformed newick tree structure"。
        if not dated_tree_newick or dated_tree_newick == ";":
            raise ResultParsingError(
                "treePL failed to produce a dated tree (empty output). "
                "This usually means the optimizer could not find feasible start "
                "rates/dates. Check that calibrations are compatible with the "
                "tree's branch lengths, or try different calibration bounds."
            )

        # 如果分支长度被修改过，需要调整回去
        if self._tree_modified and self._original_tree is not None:
            dated_tree_newick = self._adjust_dated_tree_lengths(dated_tree_newick)

        # 恢复原始叶节点名并重新标注内部校准节点
        dated_tree_newick = self._restore_tree_labels(dated_tree_newick)

        # 解析节点年龄
        node_ages = self._extract_ages_from_tree(dated_tree_newick)

        # 检查是否有结果落在约束边界（treePL 的 hard min/max 可能把解推到边界）
        boundary_warnings = self._check_constraint_boundaries(node_ages)
        self._warnings.extend(boundary_warnings)

        # 收集元数据（既装 scalar 也装 per-run CV 记录，因此值类型是 Any）
        metadata: Dict[str, Any] = {
            "optimal_smooth": self._optimal_smooth,
            "seed": self.common_config.seed,
            "cv_type": "randomcv" if self.config.randomcv else "cv",
            "cv_replicates": self.config.cv_replicates,
            "prime_replicates": self.config.prime_replicates,
            "tree_modified": self._tree_modified,
        }

        if self._prime_result and self._prime_result.success:
            metadata["opt"] = self._prime_result.opt
            metadata["optad"] = self._prime_result.optad

        if self._cv_results:
            metadata["cv_results"] = [
                {"smoothing": r.smoothing, "chisq": r.chisq}
                for r in self._cv_results
                if r.success
            ]

        result = DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=0.0,
            metadata=metadata,
        )

        for warning in self._warnings:
            result.add_warning(warning)

        self.logger.success(f"Parsed treePL results: {len(node_ages)} node ages")
        return result

    def _check_constraint_boundaries(
        self, node_ages: Dict[str, NodeAgeEstimate]
    ) -> List[Any]:
        """检查 treePL 结果是否落在约束边界上。

        treePL 只支持 hard min/max，使用 uniform/soft 约束时优化器可能把年龄
        推到边界。若检测到边界解，给出警告提示用户校准区间可能过窄。
        """
        warnings: List[Any] = []
        if not self._calibrations:
            return warnings

        tol_pct = 0.005  # 0.5% 相对容差
        tol_abs = 0.05  # 绝对容差（Ma）

        for cal in self._calibrations:
            estimate = node_ages.get(cal.name)
            if estimate is None:
                continue
            age = estimate.mean_age
            constraint = cal.age_constraint

            def _is_at_boundary(value: float, bound: float) -> bool:
                tolerance = max(tol_abs, abs(bound) * tol_pct)
                return abs(value - bound) <= tolerance

            at_boundary = False
            bound_info = ""
            if isinstance(constraint, UniformAgeConstraint):
                if _is_at_boundary(age, constraint.min_age):
                    at_boundary = True
                    bound_info = f"min={constraint.min_age}"
                elif _is_at_boundary(age, constraint.max_age):
                    at_boundary = True
                    bound_info = f"max={constraint.max_age}"
            elif isinstance(constraint, MaximumAgeConstraint):
                if _is_at_boundary(age, constraint.max_age):
                    at_boundary = True
                    bound_info = f"max={constraint.max_age}"
            elif isinstance(constraint, SoftLowerBoundConstraint):
                if _is_at_boundary(age, constraint.min_age):
                    at_boundary = True
                    bound_info = f"min={constraint.min_age}"

            if at_boundary:
                msg = (
                    f"treePL result for '{cal.name}' ({age:.3f} Ma) lies on "
                    f"constraint boundary ({bound_info}). This may indicate the "
                    f"calibration bound is too narrow for treePL's hard-constraint optimizer."
                )
                warnings.append(SemanticDegradationWarning(msg))
                self.logger.warning(msg, category=SemanticDegradationWarning)

        return warnings

    def _restore_tree_labels(self, newick: str) -> str:
        """恢复 treePL 输出树的原始叶节点名和内部校准节点标签。

        treePL 执行时使用的是被截断/映射后的短名，且内部节点标签被移除。
        本方法使用 prepare_inputs 阶段记录的 `_name_mapping` 和校准信息：
          1) 先将叶节点短名替换回原始全名；
          2) 再用原始全名（cal.mrca_leaf_pair）为各校准 MRCA 节点添加标签。
        顺序不可颠倒：一旦叶节点仍为短名，用全名定位 MRCA 会失败。
        """
        name_mapping = getattr(self, "_name_mapping", None)
        if not name_mapping:
            # 没有做过短名映射，只需重新标注内部节点
            return self._annotate_internal_nodes(newick, reverse_mapping=None)

        # 短名 -> 全名
        reverse_mapping = {short: full for full, short in name_mapping.items()}

        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)

            # 1. 先恢复叶节点名（必须在标注内部节点之前）
            for leaf in tree.get_leaves():
                full_name = reverse_mapping.get(leaf.name)
                if full_name:
                    leaf.name = full_name

            # 2. 再用原始全名重新标注内部校准节点
            for cal in self._calibrations:
                if not cal.name:
                    continue
                if cal.is_root_node:
                    tree.name = cal.name
                    continue
                if not cal.mrca_leaf_pair or len(cal.mrca_leaf_pair) < 2:
                    continue
                tip1, tip2 = cal.mrca_leaf_pair
                try:
                    ancestor = tree.get_common_ancestor(tip1, tip2)
                    ancestor.name = cal.name
                except Exception as e:
                    self.logger.debug(
                        f"Could not annotate calibration node '{cal.name}' in treePL output: {e}"
                    )

            # 确保根节点标签被写入 Newick（ETE3 默认会丢弃根节点名）
            root_name = None
            for cal in self._calibrations:
                if cal.is_root_node:
                    root_name = cal.name
                    break
            newick_out: str = tree.write(format=1)
            if root_name:
                newick_out = self._ensure_root_label_in_newick(newick_out, root_name)
            return newick_out
        except ImportError:
            self.logger.warning("ETE3 not available, cannot restore treePL tree labels")
            return self._annotate_internal_nodes(newick, reverse_mapping)
        except Exception as e:
            self.logger.warning(f"Failed to restore treePL tree labels: {e}")
            return newick

    def _ensure_root_label_in_newick(self, newick: str, root_name: str) -> str:
        """确保 Newick 字符串中包含根节点标签。

        ETE3 的 Tree.write(format=1) 会丢弃根节点名称，因此需要手动在末尾
        的分号前注入根节点名。
        """
        if not root_name or not newick:
            return newick
        newick = newick.strip()
        if not newick.endswith(";"):
            return newick
        # 如果根节点名已存在，则不做处理
        if re.search(
            rf"\){re.escape(root_name)}(?:[\d.eE+-]*)?;\s*$", newick
        ) or re.search(rf"\){re.escape(root_name)}:[\d.eE+-]+;\s*$", newick):
            return newick

        body = newick[:-1].rstrip()
        # 无根分支长度：(...); → (...)Root;
        if body.endswith(")"):
            return f"{body}{root_name};"
        # 有根分支长度：(...):X; → (...)Root:X;
        match = re.search(r"\)(:[\d.eE+-]+)$", body)
        if match:
            prefix = body[: match.start() + 1]  # 包含最后的 ')'
            suffix = match.group(1)  # 如 :0
            return f"{prefix}{root_name}{suffix};"
        return newick

    def _annotate_internal_nodes(
        self, newick: str, reverse_mapping: Optional[Dict[str, str]]
    ) -> str:
        """仅重新标注内部校准节点（不做叶节点名恢复时的回退）。"""
        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)
            for cal in self._calibrations:
                if not cal.name or cal.is_root_node:
                    continue
                if not cal.mrca_leaf_pair or len(cal.mrca_leaf_pair) < 2:
                    continue
                try:
                    ancestor = tree.get_common_ancestor(*cal.mrca_leaf_pair)
                    ancestor.name = cal.name
                except Exception as e:
                    # C-8：标不上名字只是"输出树少一个标签"，不值得中止整次运行，
                    # 但静默吞掉会让下游以为该校准没有对应节点，必须留痕。
                    self.logger.warning(
                        f"Could not annotate calibration node '{cal.name}' in the "
                        f"treePL output tree: {e}"
                    )
            # ete3 没有类型标注：write() 在运行期返回 str，这里把它落到一个
            # 带标注的局部变量上，避免把 Any 直接交给声明为 str 的返回值。
            rendered: str = tree.write(format=1)
            return rendered
        except Exception as e:
            self.logger.warning(
                f"Internal-node annotation skipped for treePL output tree: {e}"
            )
            return newick

    def _adjust_dated_tree_lengths(self, newick: str) -> str:
        """
        调整定年树的分支长度（如果之前被乘过）

        将分支长度除以相同的乘数，恢复原始比例
        """
        if not self._tree_modified:
            return newick

        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)

            # 将所有分支长度除以实际使用的乘数
            for node in tree.traverse():
                if node.dist is not None:
                    node.dist /= self._actual_multiplier

            rendered: str = tree.write(format=5)
            return rendered

        except ImportError:
            self.logger.warning("ETE3 not available, cannot adjust dated tree lengths")
            return newick

    def _extract_ages_from_tree(self, newick: str) -> Dict[str, NodeAgeEstimate]:
        """
        从定年树中提取节点年龄

        treePL 输出的定年树以分支长度表示时间（Ma），内部节点通常没有标签。
        这里通过 MRCA 叶节点对重新定位校准节点，并计算节点到根的距离来得到年龄。
        """
        node_ages = {}

        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)
            root = tree.get_tree_root()

            # 树根到最远叶节点的距离即根节点年龄
            root_age = max(root.get_distance(leaf) for leaf in tree.get_leaves())

            for cal in self._calibrations:
                # 根节点校准：直接使用根年龄
                if cal.is_root_node:
                    node_ages[cal.name] = NodeAgeEstimate(
                        mean_age=root_age, ci_type=CIType.NONE
                    )
                    continue

                # 确定用于定位 MRCA 的叶节点对
                if cal.mrca_leaf_pair:
                    tip1, tip2 = cal.mrca_leaf_pair
                elif len(cal.resolved_taxa) >= 2:
                    try:
                        from ..models import PhylogeneticTree

                        pt = PhylogeneticTree.from_newick(newick)
                        cal_pair = pt.get_mrca_terminals(cal.resolved_taxa)
                        tip1, tip2 = cal_pair
                    except Exception:
                        self.logger.warning(
                            f"Skipping calibration '{cal.name}' in treePL result extraction: "
                            f"no mrca_leaf_pair"
                        )
                        continue
                else:
                    self.logger.warning(
                        f"Skipping calibration '{cal.name}' in treePL result extraction: "
                        f"no mrca_leaf_pair and insufficient resolved_taxa"
                    )
                    continue

                try:
                    ancestor = tree.get_common_ancestor(tip1, tip2)
                    root_to_ancestor = root.get_distance(ancestor)
                    age = root_age - root_to_ancestor

                    node_ages[cal.name] = NodeAgeEstimate(
                        mean_age=age, ci_type=CIType.NONE
                    )
                except Exception as e:
                    self.logger.warning(f"Could not extract age for {cal.name}: {e}")

            # 补充未命名的内部节点，使报告能显示全部 7 个内部节点
            try:
                calibrated_node_ids = set()
                for cal in self._calibrations:
                    if cal.is_root_node:
                        calibrated_node_ids.add(id(root))
                        continue
                    if not cal.mrca_leaf_pair:
                        continue
                    tip1, tip2 = cal.mrca_leaf_pair
                    ancestor = tree.get_common_ancestor(tip1, tip2)
                    calibrated_node_ids.add(id(ancestor))

                idx = 1
                for node in tree.traverse("postorder"):
                    if node.is_leaf() or id(node) in calibrated_node_ids:
                        continue
                    age = root_age - root.get_distance(node)
                    if age < 0:
                        continue
                    node_ages[f"internal_node_{idx}"] = NodeAgeEstimate(
                        mean_age=age, ci_type=CIType.NONE
                    )
                    idx += 1
            except Exception as e:
                self.logger.warning(
                    f"Could not extract unnamed internal nodes for treePL: {e}"
                )

        except ImportError:
            self.logger.warning("ETE3 not available, cannot extract ages from tree")

        return node_ages


# 注册适配器
DatingMethodRegistry.register("treepl", TreePLMethod)
