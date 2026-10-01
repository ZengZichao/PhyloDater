"""
MCMCTreeMethod - MCMCTree (PAML) 定年适配器

实现 MCMCTree 定年功能
根据诊断报告（4.8 和 4.10.8）：
- 支持密码子模型 (seqtype=1)
- 处理蛋白质序列的 Hessian 矩阵计算
- 正确的 RootAge 格式 ('<value')
- Gamma-Dirichlet 先验 (rgene_gamma 4 参数)
- 终止密码子处理
"""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

if TYPE_CHECKING:
    import dendropy

from ..core import DatingMethod, DatingMethodRegistry
from ..core.exceptions import CalibrationError, ExecutionError
from ..infrastructure import ProcessRunner, get_logger
from ..infrastructure.checkpoint import MCMCTreeCheckpointManager
from ..infrastructure.configuration import (
    CommonConfig,
    MCMCTreeConfig,
    SoftwarePaths,
    ToolConfig,
)
from ..infrastructure.safe_io import safe_writer
from ..models import (
    AgeConstraint,
    CalibrationPoint,
    CIType,
    DatingResult,
    FixedAgeConstraint,
    NodeAgeEstimate,
    PhylogeneticTree,
    UniformAgeConstraint,
)


def _has_path_separator(path: str) -> bool:
    """判断路径字符串是否包含目录分隔符（即不是纯命令名）。"""
    return (
        ("/" in path)
        or (os.sep in path)
        or (os.altsep is not None and os.altsep in path)
    )


def _is_dat_file(path: Path) -> bool:
    """判断路径是否为 PAML 模型矩阵 ``.dat`` 文件（扩展名大小写不敏感）。"""
    return path.suffix.lower() == ".dat"


def _detect_paml_version(version_banner: Optional[str]) -> Optional[str]:
    """从 MCMCTree --version 输出中解析 PAML 版本号。

    典型输出："MCMCTREE in paml version 4.10.8, 15 May 2025"
    返回例如 "4.10.8" 或 "4.8"；无法解析时返回 None。
    """
    if not version_banner:
        return None
    match = re.search(
        r"paml\s+version\s+(\d+(?:\.\d+)+)", version_banner, re.IGNORECASE
    )
    if match:
        return match.group(1)
    # 兼容某些旧版本只输出 "PAML v4.8" 的情况
    match = re.search(r"PAML\s+v?\s*(\d+(?:\.\d+)+)", version_banner, re.IGNORECASE)
    if match:
        return match.group(1)
    return None


def _paml_versions_compatible(expected: str, actual: str) -> bool:
    """判断两个 PAML 版本字符串是否兼容（用于版本校验）。

    判据是 **主.次版本系列**（如 ``4.8`` 与 ``4.10``）：控制文件格式、数据类型
    支持与速率估算流程的差异都发生在系列之间，同系列内的补丁号（4.10.7 /
    4.10.8 / 4.10.10）不改变 PhyloDater 生成的 ``.ctl`` 语义。因此：

    - ``4.10.8`` vs ``4.10.10`` → 兼容（bioconda 当前发布的 ``paml`` 就是 4.10.10，
      要求补丁号逐字相等会让按文档装好的环境直接跑不起来）；
    - ``4.8`` vs ``4.10.x`` → 不兼容，仍然阻断，避免控制文件错版。
    """
    if not expected or not actual:
        return True  # 无法判断时不阻断

    def _series(version: str) -> Optional[Tuple[int, ...]]:
        parts = re.findall(r"\d+", version.strip())
        if not parts:
            return None
        nums = tuple(int(p) for p in parts[:2])
        return nums if len(nums) >= 2 else (nums[0], 0)

    expected_series = _series(expected)
    actual_series = _series(actual)
    if expected_series is None or actual_series is None:
        return True  # 解析失败时不阻断，交由上层日志提示
    return expected_series == actual_series


# MCMCTree 控制文件中 RootAge 的默认宽松上界（单位：Ga）。
#
# 语义与取值口径（审阅报告 C-9）：
# - 本值 **不是** 一个地质学年龄估计，而是「用户既未提供根校准、也未配置
#   root_age」时的兜底上界，作用是让 ``RootAge = '<value'`` 这一硬上界
#   不再对根年龄产生任何实质约束。
# - 它刻意大于 ``models.constraints.MAX_GEOLOGICAL_AGE_MA``（4600 Ma = 4.6 Ga）：
#   两个模块的口径不同是有意为之——constraints 层校验的是「用户输入的化石/
#   分子年龄是否物理可能」，而这里的 RootAge 只要求「不比任何可能答案更紧」。
#   MCMCTree 的时间超先在 RootAge 上界处截断；上界取 10 Ga 时该截断位于
#   任何生物学相关年龄之外，因此对 ``time`` 超先的数值稳定性无影响
#   （相对地，若把上界压到 4.6 Ga，接近地球年龄的根后验会被截断，
#   导致 ``lnPrior`` 在边界附近出现人为尖峰与采样退化）。
# - 该兜底意味着「根年龄完全由分子钟与 BD 先验决定」，因此使用时一定伴随
#   WARNING（见 ``_generate_control_file``），避免用户误以为有无根校准之别。
DEFAULT_ROOT_AGE_GA = 10.0


# --------------------------------------------------------------------------- #
# MCMC 收敛判定阈值（审阅报告 B-4）
# --------------------------------------------------------------------------- #
# PSRF（Gelman–Rubin ``sqrt(R-hat)``）判定阈值：只有 max(PSRF) <= 该值才算收敛。
# 旧实现用 1.2（甚至 2.0 仍判「已收敛」），会把教科书意义上明显未收敛的链
# 记为 is_converged=True。此处采用 MCMC 定年文献常用的宽松值 1.05
# （Gelman & Rubin 1992；严格实践取 1.01），并要求同时满足 ESS 阈值。
DEFAULT_PSRF_THRESHOLD = 1.05
# 最小有效样本量（ESS）：任一参数的 ESS 低于该值即判未收敛。
DEFAULT_ESS_THRESHOLD = 100.0
# 两条阈值均可通过配置项 ``psrf_threshold`` / ``ess_threshold`` 覆盖
# （见 ``MCMCTreeConfig``；未提供配置项时用上面的模块默认值）。

# --------------------------------------------------------------------------- #
# mcmc.txt 列识别（审阅报告 A-1 / A-2）
# --------------------------------------------------------------------------- #
# MCMCTree 在写入任何采样之前，先由 ``collectx()`` 写出一行 **非数值表头**：
#   Gen\tt_n5\tt_n6\t...\tmu[\tsigma2][\tkappa[\talpha]]\tlnL
# 对照上游源码：
#   - 4.10.8: src/mcmctree.c:2045-2110（collectx），1536-1541（InfiniteSites）
#   - 4.8:    src/mcmctree.c:1724-1741（collectx），1337-1342（InfiniteSites）
# 因此：
#   1) ``np.loadtxt`` 必须 skiprows=1，否则第一列 'Gen' 触发 ValueError
#      （A-1：旧实现在第一条链上就 return False，PSRF/ESS 从未被计算）；
#   2) ``Gen`` 列是迭代序号、各链逐元素相同，必须在诊断前剔除（A-2）。
#      它的链间方差为 0（PSRF 恒 ≈ 1，等于掺进一条永远达标的假参数），而它的
#      组内 ESS 极低（单调趋势的自相关时间极长），会把每次正常运行拖成
#      「ESS 不足」——旧实现因此在每次正常运行时都误报 "ESS may be low"。
# 列名按下表识别（不依赖列位置），兼容不同 PAML 版本/分支的写法差异。
_MCMC_ITERATION_COLUMN_NAMES = frozenset(
    {
        "gen",
        "generation",
        "iter",
        "iters",
        "iterno",
        "iteration",
        "sample",
    }
)

# 多链区间口径的如实描述（审阅报告 B-3）：R 条链各自 95% 区间的「外包络」
# （下界取最小、上界取最大）不是任何后验概率区间，因此不配 HPD95 标签。
_ENVELOPE_INTERVAL_DESCRIPTION = (
    "outer envelope across chains (min of chain lowers / max of chain uppers); "
    "NOT a 95% HPD"
)
_SINGLE_CHAIN_INTERVAL_DESCRIPTION = "single-chain 95% HPD as reported by MCMCTree"


class MCMCTreeMethod(DatingMethod[MCMCTreeConfig]):
    """
    MCMCTree (PAML) 定年适配器

    特性：
    - 两阶段执行：usedata=3 (Hessian) -> usedata=2 (MCMC)
    - 支持多种序列类型 (DNA, Codon, Protein)
    - 支持复杂先验分布 (Gamma, Skew Normal, Skew T)
    - 蛋白质序列手动修复工作流
    - 检查点/恢复功能
    """

    def __init__(
        self,
        config: Union[ToolConfig, MCMCTreeConfig],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
        enable_checkpoint: bool = True,
    ) -> None:
        # 入口保留双形态（历史上两种调用都存在），但在交给基类之前先归一：
        # 基类把 ``self.config`` 记为 ``MCMCTreeConfig``，适配器后续只读子配置。
        # ``isinstance`` 与原先的 ``hasattr(config, "mcmctree")`` 对这两种入参等价：
        # ``ToolConfig`` 必有 ``mcmctree``，``MCMCTreeConfig`` 必无。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            method_config = config.mcmctree
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
        self._has_root_calibration = False
        # PAML 4.10+ 弃用 RootAge，根约束是否已写入树文件
        self._root_annotation_in_tree = False
        self._seqtype = 0  # 0=DNA, 1=Codon, 2=Protein
        self._is_protein = False
        self._enable_checkpoint = enable_checkpoint
        self._checkpoint_manager: Optional[MCMCTreeCheckpointManager] = None
        # 实际探测到的 PAML 版本（由 validate_environment 填充；未探测为 None）
        self._detected_paml_version: Optional[str] = None
        # 阶段 0 速率估算结果: rgene_gamma 的 beta (均值 = rate_alpha/beta) 与估算速率
        self._calculated_beta: Optional[float] = None
        self._estimated_rate: Optional[float] = None
        # PAML dat/ 目录审计结果（由 _audit_paml_dat_files 填充，见 B-5）
        self._dat_dir: Optional[Path] = None
        self._missing_dat_files: List[str] = []
        self._required_dat_files_missing: List[str] = []
        # 最近一次 MCMC 收敛诊断的量化结果（PSRF/ESS/阈值/样本数），见 B-4
        self._convergence_diagnostics: Optional[Dict[str, Any]] = None

    @property
    def method_name(self) -> str:
        return "mcmctree"

    def _get_executable(self, name: str) -> str:
        """获取 PAML 可执行文件路径，优先使用 software_paths 自定义路径。

        多版本 PAML（4.8 与 4.10.8 等）的二进制通常同名（mcmctree/codeml/baseml）。
        推荐通过以下任一方式区分：
          1) 在 software_paths 中显式指定完整路径，例如
             mcmctree_bin=/opt/paml-4.10.8/bin/mcmctree
          2) 指定 paml_path=/opt/paml-4.10.8，本方法会自动拼接 bin/<name>
          3) 为不同版本创建重命名/符号链接（如 mcmctree48、mcmctree4108）
             并通过 *_bin 指定。
        """
        # 1) 用户显式指定了该二进制路径
        explicit = self.get_software_path(f"{name}_bin")
        if explicit:
            # 若显式值只是命令名（无路径分隔符），优先尝试用 paml_path 补全
            if not _has_path_separator(explicit):
                paml_path = self.get_software_path("paml_path")
                if paml_path:
                    candidate = Path(paml_path) / "bin" / explicit
                    if candidate.exists() and os.access(candidate, os.X_OK):
                        return str(candidate.resolve())
            return explicit

        # 2) 未显式指定，但提供了 paml_path，使用 <paml_path>/bin/<name>
        paml_path = self.get_software_path("paml_path")
        if paml_path:
            candidate = Path(paml_path) / "bin" / name
            if candidate.exists() and os.access(candidate, os.X_OK):
                return str(candidate.resolve())
            # 兼容 Windows 参考发行版（如 paml-4.8/bin/mcmctree.exe）
            candidate_exe = Path(paml_path) / "bin" / f"{name}.exe"
            if candidate_exe.exists():
                return str(candidate_exe.resolve())

        # 3) 回退到 PATH 中的命令名
        return name

    def _get_paml_bin_dir(self) -> Optional[str]:
        """获取应前置到 PATH 的 PAML 二进制目录。

        MCMCTree 在执行 Hessian/MCMC 阶段时会按 PATH 调用 baseml/codeml 等辅助
        二进制。当用户通过 ``paml_path`` 或 ``mcmctree_bin`` 指定版本时，需要把
        对应目录前置到 PATH，避免 4.8 的 mcmctree 调用到 4.10.8 的 baseml/codeml。
        """
        paml_path = self.get_software_path("paml_path")
        if paml_path:
            bin_dir = Path(paml_path) / "bin"
            if bin_dir.exists() and bin_dir.is_dir():
                return str(bin_dir.resolve())

        explicit = self.get_software_path("mcmctree_bin")
        if explicit and _has_path_separator(explicit):
            exe_path = Path(explicit)
            if exe_path.exists():
                return str(exe_path.parent.resolve())

        return None

    def _get_paml_env(self) -> Dict[str, str]:
        """返回将 PAML 二进制目录前置到 PATH 的环境变量副本。"""
        import os as _os

        env = _os.environ.copy()
        bin_dir = self._get_paml_bin_dir()
        if bin_dir:
            current_path = env.get("PATH", "")
            env["PATH"] = f"{bin_dir}{_os.pathsep}{current_path}"
            self.logger.debug(f"Prepended PAML bin dir to PATH: {bin_dir}")
        return env

    def validate_environment(self) -> bool:
        """检测 MCMCTree 是否可用，校验实际版本与期望版本一致，并审计 .dat 矩阵。"""
        from ..infrastructure.software_dependencies import SoftwareDependencyManager

        runner = ProcessRunner()
        executable = self._get_executable("mcmctree")

        available = runner.check_executable(executable)

        if available:
            version_banner = runner.get_version(executable, "--version")
            self.logger.info(
                f"MCMCTree detected: {version_banner or 'unknown version'}"
            )
            expected = getattr(self.config, "paml_version", "4.10.8")
            actual = _detect_paml_version(version_banner)
            # 记录产出结果的后端软件版本（写入 runtime_metadata.json）。
            self.software_version = version_banner or f"PAML {actual or 'unknown'}"
            # 保存实际探测到的 PAML 版本，供后续控制文件格式选择/调试使用
            self._detected_paml_version = actual
            self.logger.info(f"Detected PAML version: {actual or 'unknown'}")
            if actual and not _paml_versions_compatible(expected, actual):
                self.logger.error(
                    f"MCMCTree 版本不匹配: 期望 PAML {expected} 系列, "
                    f"但可执行文件 '{executable}' 实际为 PAML {actual}。"
                    f"4.8 与 4.10.x 的控制文件格式不同，因此不会自动降级。"
                    f"请通过 --method-args 'mcmctree_bin=/path/to/paml-{expected}/bin/mcmctree' "
                    f"或配置 paml_path=/path/to/paml-{expected} 指定正确版本，"
                    f"或改用 --method-args 'paml_version={actual}' 按实际安装的版本运行。"
                )
                return False
            if actual and actual.strip() != str(expected).strip():
                self.logger.info(
                    f"MCMCTree: configured PAML version {expected} differs from the "
                    f"detected {actual}, but both belong to the same release series; "
                    "continuing with the installed executable."
                )
            # B-5: .dat 模型矩阵审计。目录存在但控制文件依赖的矩阵缺失是可证明的
            # 安装缺陷 → 直接判为环境不可用（旧实现连 dat 目录不存在都只记 debug，
            # 用户要等几十分钟才在 codeml/mcmctree 里看到一句难懂的报错）。
            if not self._validate_paml_dat_files():
                return False
        else:
            self.logger.warning(f"MCMCTree not found: {executable}")
            install_cmd = SoftwareDependencyManager.get_install_command("mcmctree")
            self.logger.warning(f"Installation suggestions:\n{install_cmd}")

        return available

    def _validate_paml_dat_files(self) -> bool:
        """审计 PAML ``dat/`` 目录；返回 False 表示存在阻断性缺陷（B-5）。"""
        audit = self._audit_paml_dat_files()
        dat_dir: Optional[Path] = audit["dat_dir"]
        self._dat_dir = dat_dir

        if dat_dir is None:
            exe = self._get_executable("mcmctree")
            self.logger.warning(
                f"PAML 'dat/' directory not found near executable '{exe}'. "
                "蛋白质模型矩阵（wag.dat 等）需要在运行前被复制进工作目录；"
                "若序列类型为 protein，Hessian/MCMC 阶段可能以 'cannot open file' 失败。"
                "请将 paml_path 指向包含 dat/ 的 PAML 安装根目录。"
            )
            return True  # 无法证明安装有问题：不阻断，仅警告

        self._missing_dat_files = list(audit["missing"])
        self._required_dat_files_missing = list(audit["required_missing"])

        if audit["missing"]:
            self.logger.warning(
                f"PAML dat directory {dat_dir} is missing {len(audit['missing'])} "
                f"expected model matrix file(s): {', '.join(audit['missing'])}"
            )
        if audit["required_missing"]:
            self.logger.error(
                f"PAML dat directory {dat_dir} lacks "
                f"{', '.join(audit['required_missing'])}, which PhyloDater's generated "
                "control files reference explicitly (aaRatefile). "
                "MCMCTree 环境检查失败：请检查 PAML 安装是否完整。"
            )
            return False
        return True

    @staticmethod
    def _usable_mrca_pair(pair: Any) -> Optional[Tuple[str, ...]]:
        """把 ``mrca_leaf_pair`` 归一成可用于定位 MRCA 的叶名元组。

        ``None``（根校准等合法情形）、非序列对象与长度 < 2 的对都返回 ``None``，
        调用方据此跳过而不对 ``None`` 做下标访问。返回归一后的元组（而不是
        ``bool``/``TypeGuard``）有两个理由：一是 ``typing.TypeGuard`` 要 Python
        3.10+，本包的支持下限不能因为一个类型工具被抬高；二是调用方紧接着
        就要取 ``[0]`` / ``[1]``，拿到已收窄的值比拿到一个断言更直接。
        多叶 MRCA 同样合法，所以是 ``Tuple[str, ...]`` 而不是固定二元组。
        """
        if pair is None:
            return None
        try:
            length = len(pair)
        except TypeError:
            return None
        if length < 2:
            return None
        return tuple(str(name) for name in pair)

    def validate_calibrations(self, calibrations: List[CalibrationPoint]) -> List[str]:
        """MCMCTree 对根节点校准有特殊限制，并对「写不进树的校准」做前置拦截。

        根节点限制：MCMCTree 的 RootAge 只能设置为上界（`'<value'`），无法真正固定
        根节点年龄。如果用户提供 fixed root，运行结果会偏离预期（例如 500 Ma 的固定
        校准可能收敛到 450 Ma 左右）。因此在此明确拦截并提示用户修改。

        校准覆盖（审阅报告 B-6）：``models/tree.py`` 的写树路径对
        ``is_root_node=False`` 且 ``mrca_leaf_pair`` 为空的校准点是 **静默 continue**，
        对解析到同一 MRCA 节点的多条校准则 **后者覆盖前者**，两条都不告警。
        MCMCTree 因此会用比用户指定更少的校准跑完并给出自洽但偏了的结果，
        而报告里显示的仍是 YAML 的全套校准。这里在真正运行前把它变成显式错误。

        注意：``is_root_node=True`` 且 ``mrca_leaf_pair is None`` 是 **合法** 组合
        （根校准由控制文件的 RootAge 承载，本就不需要叶对；见审阅报告 B-20），
        因此下面的检查一律跳过根校准，绝不对 None 下标。
        """
        errors: List[str] = []
        root_cal = next((c for c in calibrations if c.is_root_node), None)
        if root_cal is not None:
            constraint = root_cal.age_constraint
            if isinstance(constraint, FixedAgeConstraint):
                errors.append(
                    f"MCMCTree 不支持固定根节点年龄。"
                    f"根节点校准 '{root_cal.name}' 当前为 fixed={constraint.fixed_age} Ma，"
                    f"但 MCMCTree 只能将其作为上界（RootAge='<{constraint.fixed_age}'）使用，"
                    f"导致根节点年龄无法精确固定。"
                )
                errors.append(
                    "建议修改方案（任选其一）：\n"
                    "  1) 将根节点约束改为 uniform，例如 min=499.9, max=500.1；\n"
                    "  2) 将根节点约束改为 maximum，例如 max=500.0；\n"
                    "  3) 如果必须精确固定根节点，请改用 treePL/pathd8/lsd2 等方法。"
                )

        # --- B-6：无法写进树的校准（静默丢点的两个来源）---------------------
        unlocatable: List[str] = []
        for cal in calibrations:
            if cal.is_root_node:
                continue
            label = cal.name or "<未命名校准>"
            pair = self._usable_mrca_pair(cal.mrca_leaf_pair)
            if pair is None:
                unlocatable.append(
                    f"  - '{label}': 非根校准缺少 mrca_leaf_pair，"
                    "写树时会被静默丢弃（MCMCTree 将不使用该校准点）。"
                )
            elif cal.age_constraint is None:
                unlocatable.append(
                    f"  - '{label}': 缺少 age_constraint，写树时会被静默丢弃。"
                )
        if unlocatable:
            errors.append(
                "以下校准点无法被写入 MCMCTree 树文件（继续运行会导致"
                "「用户 YAML 的校准数 ≠ 实际生效的校准数」）：\n"
                + "\n".join(unlocatable)
            )
            errors.append(
                "修复建议：在 YAML 中为该节点显式提供 mrca_pair（两叶节点对），"
                "或先运行校准解析（CalibrationResolver）填充 mrca_leaf_pair 后再定年。"
            )

        # 叶对完全相同的多条校准必然解析到同一 MRCA 节点 → tree.py 后者覆盖前者。
        # （不同叶对但同一节点的情况由 prepare_inputs 用树拓扑进一步核查。）
        by_pair: Dict[Tuple[str, str], List[str]] = {}
        for cal in calibrations:
            if cal.is_root_node:
                continue
            pair = self._usable_mrca_pair(cal.mrca_leaf_pair)
            if pair is None:
                continue
            key = (pair[0], pair[1])
            by_pair.setdefault(key, []).append(cal.name or "<未命名校准>")
        for pair, names in by_pair.items():
            if len(names) > 1:
                errors.append(
                    f"多条校准点使用同一 MRCA 叶对 {pair}（{', '.join(names)}），"
                    "写树时只有一条会被写入 MCMCTree 节点名，其余被静默覆盖。"
                    "请合并为一条约束（区间取下界最大、上界最小），或改用不同的叶对。"
                )

        return errors

    @staticmethod
    def _call_with_calibration_annotations(
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        include_root: bool = False,
        default_root_constraint: Optional[AgeConstraint] = None,
    ) -> PhylogeneticTree:
        """调用 ``tree.with_calibration_annotations`` 并保持向后兼容。

        旧版实现或测试用的 FakeTree 可能只接受 ``(calibrations)`` 一个位置参数；
        通过 ``inspect.signature`` 探测关键字支持情况，避免破坏既有调用方。
        """
        import inspect

        sig = inspect.signature(tree.with_calibration_annotations)
        params = sig.parameters
        if "include_root" in params and "default_root_constraint" in params:
            return tree.with_calibration_annotations(
                calibrations,
                include_root=include_root,
                default_root_constraint=default_root_constraint,
            )
        return tree.with_calibration_annotations(calibrations)

    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
        生成 MCMCTree 输入文件

        阶段 1: 名称标准化和序列类型检测
        """
        self._calibrations = calibrations

        # 检测序列类型
        if alignment_path and alignment_path.exists():
            self._detect_sequence_type(alignment_path)

        # 检查根节点校准
        self._has_root_calibration = any(cal.is_root_node for cal in calibrations)

        # B-6 前置审计：写树前核对「用户请求的校准」是否都能被真正写入。
        # 必须在名称压缩（下面的 name_mapping）之前执行——写树用的是全名树。
        self._audit_calibration_coverage(tree, calibrations)

        # 建立统一名称映射（ACCESSION 模式，≤30 字符，PHYLIP 安全）。
        # PAML 的 PHYLIP 格式截断到 30 字符，序列文件会用短名；树必须用同一套短名，
        # 否则 MCMCTree 报 "species XXX not found in main tree"。
        # 保存映射供 parse_results 做反向（短名→全名）查找：mcmctree 输出树使用短名，
        # 而校准点 mrca_leaf_pair 仍用全名，需靠此映射在时间树中定位 MRCA。
        name_mapping = self._build_name_mapping(tree)
        self._name_mapping = name_mapping  # Dict[全名, 短名] 或 None

        # 根据 PAML 版本选择树文件格式
        tree_file = self.work_dir / "mcmctree.tree"

        if self._is_paml_48():
            # PAML 4.8: 使用 NEXUS/UTREE 风格，约束内联在节点名中，
            # 并在树字符串开头自带 PHYLIP 头 "n_species 1"。
            tree_with_calib = tree.with_paml48_calibration_annotations(calibrations)
        elif self._is_paml_4108_or_later():
            # PAML 4.10+ 已弃用控制文件中的 RootAge，必须将根约束写入树。
            # 未提供根校准时写入宽松默认界（取约束层允许的最大地质年龄），
            # 避免 MCMCTree 以 "Only bounds for the root age are implemented" 中止。
            from ..models.constraints import MAX_GEOLOGICAL_AGE_MA

            default_root = UniformAgeConstraint(
                min_age=0.001, max_age=MAX_GEOLOGICAL_AGE_MA
            )
            # 保持对旧版/测试用 FakeTree 的兼容：仅当 with_calibration_annotations
            # 接受 include_root 关键字时才传入。
            tree_with_calib = self._call_with_calibration_annotations(
                tree,
                calibrations,
                include_root=True,
                default_root_constraint=default_root,
            )
            self._root_annotation_in_tree = True
        else:
            # PAML 4.9: 使用标准 Newick，根约束仍由 RootAge 承载。
            tree_with_calib = self._call_with_calibration_annotations(
                tree, calibrations
            )

        # B-6 出口核验：写树实现（models/tree.py）用 id(MRCA) 作键、并对缺少
        # mrca_leaf_pair 的校准直接 continue，因此「请求 N 条」与「实际写进树的
        # 条数」可能不等。这里按每条校准渲染出的约束串逐一回查，任何一条没出现
        # 就报错，绝不带着更少的校准继续跑。
        self._verify_calibrations_written(tree_with_calib, calibrations)

        # 把树里的全名替换成与序列文件一致的短名
        if name_mapping is not None:
            tree_with_calib = tree_with_calib.with_renamed_leaves(name_mapping)

        tree_with_calib.write(tree_file)

        # PAML 4.9+ / 4.10.8+ 需要单独追加 PHYLIP 头
        # （4.8 已在 with_paml48_calibration_annotations 中自带头部）
        if not self._is_paml_48():
            self._add_phylip_header_to_tree(tree_file, tree)

        # 准备序列文件（PHYLIP 格式，使用与树一致的短名）
        seq_file = self.work_dir / "mcmctree.phy"
        if alignment_path:
            self._convert_to_phylip(alignment_path, seq_file, name_mapping)

        # 把 PAML 安装目录下的 .dat 文件复制到工作目录，
        # 否则 codeml/mcmctree 在临时工作目录运行时找不到 wag.dat 等。
        self._copy_paml_dat_files()

        self._input_files = {"tree": tree_file, "seq": seq_file}

        self.logger.info("Generated MCMCTree input files")

        return self._input_files

    def _audit_calibration_coverage(
        self, tree: PhylogeneticTree, calibrations: List[CalibrationPoint]
    ) -> None:
        """写树前的校准覆盖审计（审阅报告 B-6），发现问题即报错而非静默丢点。

        覆盖两类静默失败：

        1. 非根校准缺少 ``mrca_leaf_pair`` —— 写树路径（``models/tree.py``）会
           直接 ``continue``，MCMCTree 于是用更少的校准跑完，输出自洽但偏了的年龄。
        2. 两条校准解析到 **同一个 MRCA 节点** —— 写树路径用 ``id(mrca)`` 作字典键，
           后者覆盖前者。用「MRCA 所包含的叶集合」判定同一节点（同一节点 ⟺ 同一
           叶集合），从而能抓住叶对不同但节点相同的嵌套分类群情形。

        拓扑判定通过 :class:`PhylogeneticTree` 的既有接口（``get_mrca``）完成，本方法
        不直接 import ete3；``get_mrca`` 在后端不可用时返回 ``None``，此时退化为
        「仅按叶对是否相同判定」，并把无法做拓扑核查这件事本身告诉用户。
        """
        structural_errors = self.validate_calibrations(calibrations)
        if structural_errors:
            # prepare_inputs 可能被直接调用（不经流水线的 validate_calibrations 闸门），
            # 因此这里的结论必须是致命的而不是只返回清单。
            message = "MCMCTree 校准配置检查未通过：\n" + "\n".join(
                f"  - {err}" for err in structural_errors
            )
            self.logger.error(message)
            raise CalibrationError(message)

        # --- 同一 MRCA 节点被多条校准占用（叶对不同但节点相同）---------------
        signature_map: Dict[Any, List[str]] = {}
        topology_check_available = False
        for cal in calibrations:
            if cal.is_root_node:
                continue
            pair = tuple(str(t) for t in (cal.mrca_leaf_pair or ())[:2])
            try:
                leaves = tree.get_mrca(list(pair))
            except Exception as e:  # 后端异常不应让整条流水线崩在审计里
                self.logger.debug(
                    f"MCMCTree calibration audit: get_mrca({pair}) failed: {e}"
                )
                leaves = None
            if leaves:
                topology_check_available = True
                signature = frozenset(leaves)
            else:
                # 无拓扑信息：用叶对本身作键（最坏情况下退化为原有的粗判）
                signature = None
            if signature is not None:
                signature_map.setdefault(signature, []).append(
                    f"{cal.name or '<未命名校准>'} (mrca_pair={pair})"
                )

        collisions = [
            (sig, names) for sig, names in signature_map.items() if len(names) > 1
        ]
        if collisions:
            detail = "\n".join(
                f"  - 节点（叶集合 {sorted(sig)}）被多条校准占用：{'; '.join(names)}"
                for sig, names in collisions
            )
            message = (
                "MCMCTree 校准冲突：多条校准解析到同一个内部节点。写树时只有最后一条"
                "会被写进节点名，其余会被静默覆盖，从而让实际生效的校准少于用户指定。\n"
                + detail
                + "\n请把这些校准合并为一条（区间下界取各下界的最大值、上界取各上界"
                "的最小值；若交集为空说明校准互相矛盾），或为不同节点选择不同的叶对。"
            )
            self.logger.error(message)
            raise CalibrationError(message)

        if not topology_check_available and len(calibrations) > 1:
            self.logger.warning(
                "MCMCTree 校准审计无法进行拓扑级同节点冲突检测"
                "（树后端不可用或 get_mrca 返回空），仅比对了叶对是否相同。"
                "若存在嵌套分类群共享同一 MRCA 节点的校准，仍可能互相覆盖。"
            )

    def _verify_calibrations_written(
        self,
        tree_with_calib: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
    ) -> None:
        """核验每条非根校准的约束串确实出现在写出的树字符串里（B-6 出口检查）。

        这是对写树实现（含其任何未来分支）的**黑盒对账**：即便上游用别的方式
        丢了点，这里也会发现。根校准不写进树（由控制文件的 ``RootAge`` 或 4.8 通路
        的 ``'<value'`` 后缀承载），因此不参与本项计数。
        """
        newick = getattr(tree_with_calib, "newick", None)
        if not isinstance(newick, str) or not newick.strip():
            # 拿不到树文本（例如被 mock 替换）时不做判定，避免误报
            self.logger.debug(
                "MCMCTree calibration write verification skipped: no Newick text available."
            )
            return

        expected: Dict[str, int] = {}
        labels: Dict[str, List[str]] = {}
        unverifiable: List[str] = []
        for cal in calibrations:
            if cal.is_root_node:
                continue
            constraint = cal.age_constraint
            if constraint is None:
                continue
            try:
                rendered = constraint.to_mcmctree_calib_string()
            except Exception as e:
                rendered = ""
                self.logger.debug(
                    f"MCMCTree calibration verification: cannot render {cal.name}: {e}"
                )
            if not rendered:
                unverifiable.append(cal.name or "<未命名校准>")
                continue
            expected[rendered] = expected.get(rendered, 0) + 1
            labels.setdefault(rendered, []).append(cal.name or "<未命名校准>")

        dropped = [s for s, n in expected.items() if newick.count(s) < n]
        if dropped:
            detail = "\n".join(
                f"  - {', '.join(labels[s])} → 约束 {s}（树中出现 {newick.count(s)} 次，"
                f"期望 {expected[s]} 次）"
                for s in dropped
            )
            message = (
                "MCMCTree 输入树缺失用户指定的校准标注：以下约束未被写入 "
                "mcmctree.tree（写树时被静默丢弃或同节点被覆盖）。\n"
                + detail
                + "\n拒绝在「实际校准数 < 请求校准数」的情况下继续运行。"
            )
            self.logger.error(message)
            raise CalibrationError(message)

        if unverifiable:
            self.logger.warning(
                "MCMCTree: 无法核验以下校准点是否写入树（其约束未渲染出 MCMCTree 校准串）："
                f"{', '.join(unverifiable)}。请人工检查 mcmctree.tree 的节点标注。"
            )

    def _build_name_mapping(self, tree: PhylogeneticTree) -> Optional[Dict[str, str]]:
        """为 PHYLIP/PAML 建立全名 → 安全短名的映射。

        PAML 的 PHYLIP 格式把序列名截断到 30 字符，树文件必须使用同一套短名。
        使用 NameMappingManager (ACCESSION 模式) 提取 ``GB_GCA_<accession>`` 部分，
        既短又唯一。如果所有名本身已 ≤30 字符且 PHYLIP 安全，则无需映射。
        """
        from ..services.name_mapping import NameMappingManager, NameShortenMode

        mapping = NameMappingManager()
        any_changed = False
        for original in tree.tip_names:
            short = mapping.add_entry(original, method=NameShortenMode.ACCESSION)
            if short != original:
                any_changed = True
        return mapping.get_mapping_dict() if any_changed else None

    # PAML 程序（codeml/mcmctree）在工作目录运行，需要 .dat 文件（如 wag.dat）
    # 在当前目录或 PATH 能找到。这里把它们从安装目录复制过来。
    #
    # 本清单是「与上游逐字核对过的期望文件名」，不是复制依据——复制走
    # ``_copy_paml_dat_files`` 的整目录枚举，清单只用于对账（见审阅报告 B-5）。
    # 名单来源：随附 PAML 4.8 与 4.10.8 的 ``dat/`` 目录，两版本内容完全一致
    # （18 个文件），文件名大小写以下游实际文件为准：
    #   - ``mtArt.dat``（旧清单误写 ``mtART.dat``）
    #   - ``MtZoa.dat``（旧清单误写 ``mtZOA.dat``）
    #   - 补入旧清单遗漏的 ``lg.dat``（LG 为常用蛋白模型）与 ``miyata.dat``
    #   - 删除上游根本不存在的 ``mtMet.dat`` / ``cpREV45.dat`` /
    #     ``aamodel.dat`` / ``figual.dat``
    # macOS 的 APFS 默认大小写不敏感，写错的条目标在本地「看起来正常」，
    # 而在 Linux 上会静默复制不到 → 几十分钟后 codeml/mcmctree 报
    # "cannot open file"。因此这里同时要求按真实大小写命名。
    _PAML_DAT_FILES = [
        "wag.dat",
        "lg.dat",
        "jones.dat",
        "jones-dcmut.dat",
        "dayhoff.dat",
        "dayhoff-dcmut.dat",
        "miyata.dat",
        "grantham.dat",
        "g1974a.dat",
        "g1974c.dat",
        "g1974p.dat",
        "g1974v.dat",
        "mtREV24.dat",
        "mtmam.dat",
        "mtArt.dat",
        "MtZoa.dat",
        "cpREV64.dat",
        "cpREV10.dat",
    ]

    # 由本适配器生成的控制文件 **显式引用** 的 .dat 文件（缺失即无法运行）。
    # 见 ``_generate_control_file`` / ``_generate_rate_estimation_ctl`` 中的
    # ``aaRatefile = wag.dat``，以及 ``_manual_protein_fix`` 补写的同名参数。
    _PAML_REQUIRED_DAT_FILES = ["wag.dat"]

    def _find_paml_dat_dir(self) -> Optional[Path]:
        """定位 PAML 安装目录下的 ``dat/`` 目录；找不到时返回 None。

        从 mcmctree 可执行文件路径反推：``bin/mcmctree`` → ``../dat``（PAML 官方
        发行结构），并兼容 conda 的 ``share/paml/dat`` 布局。
        """
        exe = self._get_executable("mcmctree")
        # _get_executable 可能返回 basename（如 "mcmctree"），需要解析为完整路径。
        # 用 shutil.which 在 PATH 中查找，避免 Path.resolve() 把相对路径解析为 cwd 相对。
        resolved = shutil.which(exe)
        exe_path = Path(resolved).resolve() if resolved else Path(exe).resolve()

        candidates = [
            exe_path.parent.parent / "dat",
            exe_path.parent.parent / "share" / "paml" / "dat",
            exe_path.parent / "dat",
        ]
        for candidate in candidates:
            try:
                if candidate.is_dir() and any(
                    _is_dat_file(p) for p in candidate.iterdir()
                ):
                    return candidate
            except OSError as e:  # 权限/竞态：不让探测本身炸掉流程
                self.logger.debug(
                    f"Cannot inspect candidate PAML dat dir {candidate}: {e}"
                )
        return None

    def _audit_paml_dat_files(self) -> Dict[str, Any]:
        """把 PAML ``dat/`` 目录内容与期望清单做**大小写不敏感**对账（B-5）。

        Returns:
            ``{"dat_dir": Path | None, "actual": {小写名: 真实文件名},
               "missing": [期望但实际不存在的项], "case_fixed": {清单名: 真实名},
               "required_missing": [控制文件真正依赖但缺失的项]}``

        注意：``dat_dir`` 定位失败时 ``missing`` / ``required_missing`` 保持为空——
        「找不到目录」与「目录里少文件」是两类问题，前者由调用方按 warning 处理
        （可能只是 mcmctree 不在标准安装结构里），后者才是可证明的安装缺陷。
        """
        dat_dir = self._find_paml_dat_dir()
        report: Dict[str, Any] = {
            "dat_dir": dat_dir,
            "actual": {},
            "missing": [],
            "case_fixed": {},
            "required_missing": [],
        }
        if dat_dir is None:
            return report

        try:
            actual_entries = sorted(
                p.name for p in dat_dir.iterdir() if _is_dat_file(p) and p.is_file()
            )
        except OSError as e:
            self.logger.warning(f"Cannot read PAML dat directory {dat_dir}: {e}")
            return report

        actual = {name.lower(): name for name in actual_entries}
        report["actual"] = actual

        for expected in self._PAML_DAT_FILES:
            real_name = actual.get(expected.lower())
            if real_name is None:
                report["missing"].append(expected)
            elif real_name != expected:
                # 清单与磁盘大小写不一致（Linux 上会直接找不到）
                report["case_fixed"][expected] = real_name

        for required in self._PAML_REQUIRED_DAT_FILES:
            if required.lower() not in actual:
                report["required_missing"].append(required)

        return report

    def _copy_paml_dat_files(self) -> None:
        """把 PAML 安装目录下的 **全部** .dat 文件复制到工作目录（B-5）。

        codeml/mcmctree 在临时工作目录运行，不会自动找到 conda 环境
        或 PAML 安装目录里的 .dat 文件。需要显式复制。

        与旧实现的三点差异：
        1. 不再按硬编码清单逐个 ``exists()`` 探测（缺失即静默跳过），而是**整目录
           枚举复制**，保留源文件真实大小写，因此 Linux/macOS 行为一致；
        2. 期望清单与磁盘内容做大小写不敏感对账，缺项记 WARNING、控制文件真正
           依赖的项缺失记 ERROR（旧实现连 dat 目录找不到都只记 debug 日志）；
        3. 结果写入 ``self._missing_dat_files`` 等属性，供 ``validate_environment``
           与报告层复述。
        """
        audit = self._audit_paml_dat_files()
        dat_dir: Optional[Path] = audit["dat_dir"]
        self._dat_dir = dat_dir
        self._missing_dat_files = list(audit["missing"])
        self._required_dat_files_missing = list(audit["required_missing"])

        if dat_dir is None:
            exe = self._get_executable("mcmctree")
            self.logger.warning(
                f"PAML 'dat/' directory not found near executable '{exe}'. "
                f"Protein replacement-matrix files ({', '.join(self._PAML_REQUIRED_DAT_FILES)}) "
                "will NOT be available in the working directory; codeml/mcmctree runs with "
                "seqtype=2 (or the protein manual-fix workflow) will fail with an opaque "
                "'cannot open wag.dat' style error. Install PAML with its dat/ directory, "
                "or point paml_path=/path/to/paml-x.y.z at the installation root."
            )
            return

        for mismatched, real_name in sorted(audit["case_fixed"].items()):
            self.logger.debug(
                f"PAML dat filename case differs from the expected list: "
                f"'{mismatched}' -> using real name '{real_name}'"
            )

        copied = 0
        failed: List[str] = []
        for real_name in sorted(audit["actual"].values()):
            src = dat_dir / real_name
            dst = self.work_dir / real_name
            if dst.exists():
                continue
            try:
                shutil.copy2(str(src), str(dst))
                copied += 1
            except Exception as e:
                failed.append(real_name)
                self.logger.warning(f"Failed to copy {real_name} from {dat_dir}: {e}")

        if failed:
            self.logger.error(
                f"Failed to copy {len(failed)} PAML dat file(s) into {self.work_dir}: "
                f"{', '.join(failed)}"
            )

        if self._missing_dat_files:
            self.logger.warning(
                f"PAML dat directory {dat_dir} does not contain {len(self._missing_dat_files)} "
                f"expected model matrix file(s): {', '.join(self._missing_dat_files)}. "
                "Analyses that reference them (e.g. aaRatefile) will fail."
            )
        if self._required_dat_files_missing:
            self.logger.error(
                f"PAML control files generated by PhyloDater reference "
                f"{', '.join(self._required_dat_files_missing)}, but it is missing from "
                f"{dat_dir}. The MCMCTree/codeml run will fail. Check the PAML installation "
                "or the version pointed to by paml_path/mcmctree_bin."
            )
        self.logger.info(
            f"Copied {copied} PAML .dat model matrix file(s) from {dat_dir} "
            f"into {self.work_dir}"
        )

    def _add_phylip_header_to_tree(
        self, tree_file: Path, tree: PhylogeneticTree
    ) -> None:
        """
        为树文件添加 PHYLIP 头

        根据 PAML 文档，树文件应该有 "n_species n_trees" 的头部
        例如: "4 1" 表示 4 个物种，1 棵树

        注意: 使用原始树的叶节点数量，因为 MCMCTree 的树文件格式
        不是标准 Newick，带校准标注的树会被错误解析
        """
        # 使用原始树的叶节点数量
        n_species = tree.num_tips
        n_trees = 1  # 目前只支持单棵树

        # 读取现有树内容
        with open(tree_file, "r", encoding="utf-8") as f:
            tree_content = f.read().strip()

        # 检查是否已有 PHYLIP 头
        lines = tree_content.split("\n")
        if len(lines) > 0:
            first_line = lines[0].strip()
            # 检查第一行是否符合 "n m" 格式（两个数字）
            parts = first_line.split()
            if len(parts) == 2 and all(p.isdigit() for p in parts):
                # 已有 PHYLIP 头，不需要修改
                return

        # 添加 PHYLIP 头
        with safe_writer(tree_file, encoding="utf-8", newline="\n") as f:
            f.write(f"{n_species} {n_trees}\n")
            f.write(tree_content)
            f.write("\n")

    def _detect_sequence_type(self, alignment_path: Path) -> None:
        """检测序列类型

        注意: MCMCtree 不支持密码子数据 (seqtype=1)，仅支持核苷酸 (seqtype=0) 和氨基酸 (seqtype=2)
        """
        from ..infrastructure import AlignmentMetadataExtractor

        extractor = AlignmentMetadataExtractor()

        # 检查是否为 PAML 4.10.8+，该版本不支持密码子数据
        is_paml_4108_or_later = self._is_paml_4108_or_later()

        # 自动检测
        if self.config.seqtype == "auto":
            seq_type = extractor.validate_sequence_type(alignment_path)
            if seq_type == "protein":
                self._seqtype = 2
                self._is_protein = True
            else:
                self._seqtype = 0  # DNA
        elif self.config.seqtype == "codon":
            # MCMCtree 4.10.8+ 不支持密码子数据
            if is_paml_4108_or_later:
                raise CalibrationError(
                    "MCMCtree 4.10.8+ does not support codon data (seqtype=1). "
                    "Please use nucleotide (seqtype=0) or amino acid (seqtype=2) data. "
                    "For codon-based analysis, please use CODEML or BASEML directly."
                )
            self._seqtype = 1
            # 验证序列长度是 3 的倍数
            metadata = extractor.extract(alignment_path)
            if metadata.sequence_length % 3 != 0:
                raise CalibrationError(
                    f"Sequence length ({metadata.sequence_length}) is not divisible by 3 for codon model"
                )
        elif self.config.seqtype == "protein":
            self._seqtype = 2
            self._is_protein = True
        else:
            self._seqtype = 0  # DNA

        self.logger.info(
            f"Detected sequence type: {self.config.seqtype} (seqtype={self._seqtype})"
        )

    def _convert_to_phylip(
        self,
        alignment_path: Path,
        output_path: Path,
        name_mapping: Optional[Dict[str, str]] = None,
    ) -> None:
        """
        转换为 PHYLIP 格式

        支持多分区比对 (Partitioned Alignments)，当 ndata > 1 时，
        在序列文件首行使用 G 选项，并在第二行指明各分区的长度信息。
        使用生成器模式避免一次性加载所有序列到内存。

        遵循 PAML 文档规范:
        - 序列名称最多 30 字符，不包含特殊符号
        - 使用两个连续空格作为名称结束标志
        - 终止密码子 (*) 替换为缺口 (-) 保持序列长度一致

        Args:
            name_mapping: 全名 → 短名映射（与树文件一致）。如果提供则优先使用，
                         否则回退到正则截断（≤30 字符，特殊字符替换为下划线）。
        """
        from Bio import SeqIO

        sequences = SeqIO.parse(alignment_path, "fasta")
        first_record = next(sequences)
        seq_len = len(first_record.seq)

        # 将生成器转换为列表后添加首条记录
        all_sequences = [first_record] + list(sequences)
        n_seqs = len(all_sequences)

        with safe_writer(output_path, encoding="utf-8", newline="\n") as f:
            # 处理多分区格式
            if self.config.ndata > 1:
                # 多分区格式（PAML 规范）:
                #   第一行: n_seqs <总列数> G
                #   第二行: G {ndata} <各分区列数...>
                # 说明：当前输入为单条拼接比对，按等分切分为 ndata 个分区，
                # 余数归入最后一个分区，保证各分区长度之和 == 总列数。
                # 若需真实不等长分区，需由上游传入分区元数据（超出本次最小改动范围）。
                base_len = seq_len // self.config.ndata
                remainder = seq_len % self.config.ndata

                partition_lengths = []
                for i in range(self.config.ndata):
                    # 最后一个分区包含余数
                    if i == self.config.ndata - 1:
                        partition_lengths.append(base_len + remainder)
                    else:
                        partition_lengths.append(base_len)

                # 第一行: n_seqs seq_len G
                f.write(f"{n_seqs} {seq_len} G\n")
                # 第二行: G ndata len1 len2 ...
                lengths_str = " ".join(str(length) for length in partition_lengths)
                f.write(f"G {self.config.ndata} {lengths_str}\n")
            else:
                # 单分区格式
                f.write(f"{n_seqs} {seq_len}\n")

            # 写入序列
            for record in all_sequences:
                original_id = record.id
                if name_mapping and original_id in name_mapping:
                    # 使用与树文件一致的短名（已是 PHYLIP 安全的 ≤30 字符）
                    name = name_mapping[original_id]
                else:
                    # 回退: 安全化名称并截断到 30 字符
                    name = re.sub(r"[^A-Za-z0-9_]", "_", original_id)[:30]
                # 确保名称中不包含连续空格（这会干扰 PAML 的名称解析）
                name = name.replace("  ", "_")
                # 处理终止密码子: 将 * 替换为 - (缺口)
                # 直接删除会导致序列长度不一致，违反 PHYLIP 格式
                seq = str(record.seq).replace("*", "-")
                # 使用双空格作为名称结束标志（PAML 文档规范）
                # 名称左对齐，不足 30 字符用空格填充，然后双空格分隔
                f.write(f"{name:<30}  {seq}\n")

    def execute(self, resume: bool = False) -> bool:
        """
        执行 MCMCTree 分析 (PAML 4.10.8 近似似然法)

        两阶段工作流 (PAML 4.10.8+):
        阶段 1: Hessian 矩阵计算 (usedata=3)
          - MCMCTree 调用 BASEML (核苷酸) 或 CODEML (氨基酸)
          - 生成 out.BV 文件 (分支长度, 梯度, Hessian)
        阶段 2: MCMC 采样 (usedata=2)
          - 使用 in.BV 文件进行近似似然计算
          - 执行 MCMC 采样估计分歧时间

        注意: 阶段 0 (速率估算) 默认对所有 PAML 版本执行；
              可通过 skip_rate_estimation=True 跳过

        Args:
            resume: 是否从检查点恢复
        """
        # 初始化检查点管理器
        if self._enable_checkpoint:
            self._checkpoint_manager = MCMCTreeCheckpointManager(self.work_dir)

            if resume:
                progress = self._checkpoint_manager.get_progress()
                self.logger.info(
                    f"Resuming from checkpoint: {progress['completed_steps']}/{progress['total_steps']} steps completed"
                )

        # 阶段 0: 速率估算 (默认执行; 可通过 skip_rate_estimation=True 跳过)
        # 根据教程: "We will estimate overall substitution rate in this dataset.
        # This estimates will be used in rgene_gamma setting in mcmctree.ctl later."
        # PAML 4.8 与 4.10.8 均支持此步骤 (4.9i 缺失该输出属 bug, 4.9j 已修复)
        if not self._should_skip_rate_estimation() and self._should_run_step(
            "rate_estimation"
        ):
            self.logger.step(0, "Rate estimation (baseml/codeml, clock=1)")
            self._run_with_checkpoint(
                "rate_estimation", self._run_rate_estimation_stage
            )
        elif not self._should_skip_rate_estimation():
            self.logger.info("Skipping rate estimation (already completed)")
        else:
            self.logger.info(
                "Skipping rate estimation (skip_rate_estimation=True): "
                "using default rgene_gamma prior"
            )

        # 阶段 1: Hessian 矩阵计算 (usedata=3)
        if self._should_run_step("hessian_calculation"):
            self.logger.step(1, "Hessian matrix calculation (usedata=3)")
            self._run_with_checkpoint("hessian_calculation", self._run_hessian_stage)
        else:
            self.logger.info("Skipping Hessian calculation (already completed)")

        # 阶段 2: MCMC 采样 (usedata=2)
        if self._should_run_step("mcmc_sampling"):
            self.logger.step(2, "MCMC sampling (usedata=2)")
            self._run_with_checkpoint("mcmc_sampling", self._run_mcmc_stage)
        else:
            self.logger.info("Skipping MCMC sampling (already completed)")

        return True

    def _should_run_step(self, step_name: str) -> bool:
        """检查是否需要执行步骤"""
        if not self._enable_checkpoint or not self._checkpoint_manager:
            return True
        return not self._checkpoint_manager.is_step_completed(step_name)

    def _run_with_checkpoint(self, step_name: str, func: Callable[[], Any]) -> Any:
        """带检查点的执行"""
        if not self._enable_checkpoint or not self._checkpoint_manager:
            return func()

        return self._checkpoint_manager.execute_with_checkpoint(
            step_name, func, skip_if_completed=True
        )

    def _run_rate_estimation_stage(self) -> None:
        """阶段 0: 速率估算 (baseml/codeml, clock=1)

        按 MCMCTree 教程的"粗略速率估算"步骤执行：
        1. 从带校准的 MCMCTree 树生成速率估算树（去除校准标注，
           根部追加绝对年龄标注 @age，无 @ 标注时 baseml/codeml 不输出速率）；
        2. DNA 用 baseml、蛋白质用 codeml，以全局钟 (clock=1) 拟合数据；
        3. 从输出文件提取 "Substitution rate is per time unit" 的总体替换速率；
        4. 以估算速率作为 rgene_gamma 先验均值：
           beta = rate_alpha / (rate * rate_scale)，即 Gamma 均值 = 调整后速率。

        版本行为说明（详见使用文档）：
        - PAML 4.8 与 4.10.8 均输出该速率行；4.9i 因 bug 缺失该输出，4.9j 已修复。
        - 4.10.8 在速率标题与数值之间多打印一个空行，提取正则已兼容。
        - 树中无绝对年龄标注 (@) 时两个版本都不会输出速率行。

        失败回退：无法确定根年龄或提取速率失败时，回退为默认弱信息先验
        rgene_gamma 均值 = rate_alpha/20（默认 G(2, 20)），并记录 WARNING。
        """
        rate_dir = self.work_dir / "rate_estimation"
        rate_dir.mkdir(exist_ok=True)

        # 把 PAML .dat 文件复制到速率估算子目录（codeml 蛋白模型需要 wag.dat 等）
        for dat in self.work_dir.glob("*.dat"):
            shutil.copy(dat, rate_dir / dat.name)

        # 1) 确定根年龄（Ga，与 MCMCTree 树内校准同单位）
        root_age_ga = self._resolve_root_age_ga_for_rate_estimation()
        if root_age_ga is None:
            self.logger.warning(
                "No absolute root age available for rate estimation "
                "(no root calibration and no root_age configured). "
                "Falling back to default rgene_gamma prior."
            )
            self._set_default_rate_prior()
            return

        # 2) 准备速率估算树（去除校准标注，根部加 @age）
        seq_file = self._input_files["seq"]
        rate_tree = rate_dir / "rate_est.tree"
        if not self._prepare_rate_estimation_tree(
            self._input_files["tree"], rate_tree, root_age_ga
        ):
            self._set_default_rate_prior()
            return

        shutil.copy(seq_file, rate_dir / Path(seq_file).name)

        # 3) 生成控制文件并运行 (DNA -> baseml, 蛋白 -> codeml)
        exe = "codeml" if self._seqtype == 2 else "baseml"
        self._generate_rate_estimation_ctl(rate_dir, exe, root_age_ga)

        self.logger.info(
            f"Running {exe} (clock=1) for rate estimation "
            f"(root age @{root_age_ga:.4g} Ga)..."
        )
        runner = ProcessRunner(
            cwd=rate_dir,
            timeout=self.common_config.timeout,
            env=self._get_paml_env(),
        )
        result = runner.run([self._get_executable(exe), f"{exe}.ctl"])
        if result.returncode != 0:
            self.logger.warning(
                f"{exe} exited with code {result.returncode}: {result.stderr[:200]}. "
                "Falling back to default rgene_gamma prior."
            )
            self._set_default_rate_prior()
            return

        # 4) 提取总体替换速率
        outfile = rate_dir / ("mlc" if self._seqtype == 2 else "mlb")
        rate = self._extract_substitution_rate(outfile)
        if rate is None:
            self.logger.warning(
                f"Could not extract substitution rate from {outfile.name}. "
                "Falling back to default rgene_gamma prior."
            )
            self._set_default_rate_prior()
            return

        scaled_rate = rate * self.config.rate_scale
        self._calculated_beta = self.config.rate_alpha / scaled_rate
        self._estimated_rate = rate
        self.logger.info(
            f"Estimated overall substitution rate: {rate:.6g} per time unit "
            f"(time unit = Ga, consistent with tree calibrations). "
            f"rgene_gamma = G({self.config.rate_alpha}, {self._calculated_beta:.4g}), "
            f"prior mean = {scaled_rate:.6g}"
        )

        # 保存速率供用户查阅（工作目录运行后会被清理，复制一份到输出目录留存）
        rate_summary = (
            f"substitution_rate_per_time_unit = {rate:.8g}\n"
            f"time_unit = Ga\n"
            f"root_age_ga = {root_age_ga:.4g}\n"
            f"rgene_gamma = {self.config.rate_alpha} {self._calculated_beta:.6g}\n"
        )
        (rate_dir / "estimated_rate.txt").write_text(rate_summary, encoding="utf-8")
        try:
            (self.output_dir / "mcmctree_rate_estimation.txt").write_text(
                rate_summary, encoding="utf-8"
            )
        except OSError as e:
            self.logger.warning(f"Could not save rate estimation summary: {e}")

    def _set_default_rate_prior(self) -> None:
        """回退为默认弱信息先验 rgene_gamma（均值 = rate_alpha/20，默认 G(2, 20)）"""
        self._calculated_beta = self.config.rate_alpha / 20.0
        self.logger.info(
            "Using default rgene_gamma = "
            f"G({self.config.rate_alpha}, {self._calculated_beta})"
        )

    def _resolve_root_age_ga_for_rate_estimation(self) -> Optional[float]:
        """为速率估算确定根年龄，返回 Ga 单位的浮点值；无法确定时返回 None。

        优先级：根节点校准约束的代表性年龄 > 配置项 root_age (Ma)。
        MCMCTree 树内校准以 Ga 写入（Ma ÷ 1000），因此 @ 标注必须同为 Ga，
        这样 baseml/codeml 估算出的速率才与 MCMCTree 的时间单位一致。
        """
        for cal in self._calibrations:
            if not cal.is_root_node or cal.age_constraint is None:
                continue
            age_ma = self._representative_age_ma(cal.age_constraint)
            if age_ma and age_ma > 0:
                return age_ma / 1000.0

        if self.config.root_age and self.config.root_age > 0:
            return self.config.root_age / 1000.0

        return None

    @staticmethod
    def _representative_age_ma(constraint: Optional[AgeConstraint]) -> Optional[float]:
        """从根校准约束中取代表性年龄（Ma），供速率估算的 @ 标注使用"""
        fixed_age = getattr(constraint, "fixed_age", None)
        if fixed_age:
            return float(fixed_age)

        min_age = getattr(constraint, "min_age", None)
        max_age = getattr(constraint, "max_age", None)
        if min_age and max_age:
            return (float(min_age) + float(max_age)) / 2.0
        if max_age:
            return float(max_age)
        if min_age:
            return float(min_age)

        # Gamma/Skew 等连续先验：用 alpha/beta 的均值或 location 作为代表值
        alpha = getattr(constraint, "alpha", None)
        beta = getattr(constraint, "beta", None)
        offset = getattr(constraint, "offset", 0.0) or 0.0
        if (
            isinstance(alpha, (int, float))
            and isinstance(beta, (int, float))
            and not isinstance(alpha, bool)
            and not isinstance(beta, bool)
            and beta > 0
            and alpha / beta > 1
        ):
            # MCMCTree Gamma 先验的 beta 单位为 Ga，均值 alpha/beta (Ga) → 转 Ma
            return float(alpha) / float(beta) * 1000.0 + float(offset)
        location = getattr(constraint, "location", None)
        if location:
            return float(location)
        return None

    def _prepare_rate_estimation_tree(
        self, src_tree: Path, dest_tree: Path, root_age_ga: float
    ) -> bool:
        """生成 baseml/codeml 速率估算用树。

        输入树是 MCMCTree 校准树（无分支长度，带引号包裹的校准标注）。
        处理：保留 "n_species n_trees" 头部（baseml/codeml 必需，缺失时报
        "error: end of tree file"），去除引号校准标注，
        在根部追加绝对年龄标注 @age（Ga）。clock=1 下 baseml/codeml
        以该年龄把相对时间换算为绝对时间，并输出每时间单位的替换速率。
        """
        content = src_tree.read_text(encoding="utf-8")
        lines = content.strip().split("\n")

        # 分离 "n_species n_trees" 头部与树体
        header = None
        if lines and re.match(r"^\d+\s+\d+\s*$", lines[0].strip()):
            header = lines[0].strip()
            tree_body = "\n".join(lines[1:]).strip()
        else:
            tree_body = content.strip()

        # 去除引号包裹的校准标注（如 'B(0.004, 0.006, 0.01)' / '<0.1'）
        tree_body = re.sub(r"'[^']*'", "", tree_body).strip()

        if not tree_body.startswith("("):
            self.logger.warning(
                "Rate estimation tree could not be parsed; "
                "falling back to default rgene_gamma prior."
            )
            return False

        tree_body = tree_body.rstrip(";").rstrip() + f" @{root_age_ga:.6g};"
        out = (header + "\n" if header else "") + tree_body + "\n"
        dest_tree.write_text(out, encoding="utf-8")
        self.logger.debug(f"Rate estimation tree saved: {dest_tree}")
        return True

    def _generate_rate_estimation_ctl(
        self, rate_dir: Path, exe: str, root_age_ga: float
    ) -> None:
        """生成 baseml/codeml 速率估算控制文件（clock=1，全局钟粗估）"""
        # 控制文件里的 model 既可能是用户给的字符串（直接写进去），也可能是
        # "auto" 解析出的 PAML 序号（0/7=核酸、2=蛋白质），所以是 str|int。
        model_value: Union[str, int] = self.config.model
        if model_value == "auto":
            model_value = 2 if self._seqtype == 2 else 7

        if self._seqtype == 2:
            ctl_content = f"""      seqfile = {Path(self._input_files["seq"]).name}
     treefile = rate_est.tree
      outfile = mlc
        noisy = 0
      verbose = 1
      runmode = 0
      seqtype = 2
        model = {model_value}
  aaRatefile = wag.dat
        clock = 1
    fix_omega = 0
        omega = 0.4
    fix_alpha = 0
        alpha = 0.5
       Malpha = 0
        ncatG = 4
        getSE = 0
   Small_Diff = 5e-7
"""
        else:
            ctl_content = f"""      seqfile = {Path(self._input_files["seq"]).name}
     treefile = rate_est.tree
      outfile = mlb
        noisy = 0
      verbose = 1
      runmode = 0
        model = {model_value}
        Mgene = 0
        clock = 1
    fix_kappa = 0
        kappa = 5
    fix_alpha = 0
        alpha = 0.5
       Malpha = 0
        ncatG = 4
        nhomo = 0
        getSE = 0
 RateAncestor = 0
   Small_Diff = 7e-6
"""
        (rate_dir / f"{exe}.ctl").write_text(ctl_content, encoding="utf-8")

    _RATE_PATTERN = re.compile(
        r"Substitution rate is per time unit\s*[:\s]*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
    )

    def _extract_substitution_rate(self, outfile: Path) -> Optional[float]:
        """从 baseml (mlb) / codeml (mlc) 输出提取总体替换速率。

        PAML 4.8 的数值紧跟标题行；4.10.8 在标题与数值之间多一个空行，
        正则中的 \\s* 同时兼容两种格式。
        """
        if not outfile.exists():
            return None
        match = self._RATE_PATTERN.search(
            outfile.read_text(encoding="utf-8", errors="ignore")
        )
        if match:
            return float(match.group(1))
        return None

    def _run_hessian_stage(self) -> None:
        """运行 Hessian 矩阵计算阶段 (usedata = 3)

        PAML 4.10.8 approximate likelihood workflow:
        1. usedata = 3: MCMCTree calls BASEML (nucleotides) or CODEML (amino acids)
           to generate branch lengths, gradient, and Hessian -> out.BV
        2. Copy out.BV to in.BV for MCMC stage
        """
        ctl_file = self.work_dir / "mcmctree_hessian.ctl"
        self._generate_control_file(ctl_file, usedata=3)

        cmd = [self._get_executable("mcmctree"), ctl_file.name]
        runner = ProcessRunner(
            cwd=self.work_dir,
            timeout=self.common_config.timeout or 7200,
            env=self._get_paml_env(),
        )
        result = runner.run(cmd)

        if result.returncode != 0:
            if self._is_protein and "tmp" in result.stderr:
                self.logger.warning(
                    "MCMCTree failed for protein, attempting manual fix"
                )
                self._manual_protein_fix()
            else:
                raise ExecutionError(
                    f"Hessian stage failed (returncode={result.returncode}): {result.stderr}"
                )

        out_bv = self.work_dir / "out.BV"
        in_bv = self.work_dir / "in.BV"

        if not out_bv.exists():
            raise ExecutionError(
                "Hessian stage did not generate out.BV. This may indicate a problem "
                "with the alignment file, tree file, or PAML installation. "
                "Check the output file for errors."
            )

        shutil.copy(out_bv, in_bv)
        self.logger.info(
            f"Generated in.BV from Hessian calculation (size: {in_bv.stat().st_size} bytes)"
        )

    def _cleanup_residual_tmp_ctl_files(self) -> None:
        """
        清理工作目录下残留的 tmp_ctl 临时控制文件

        在生成新的控制文件之前，清理上次失败任务可能遗留的 tmp*.ctl 文件
        避免旧文件干扰新的分析
        """
        import glob
        import os

        tmp_ctl_pattern = str(self.work_dir / "tmp*.ctl")
        tmp_ctls = glob.glob(tmp_ctl_pattern)

        if tmp_ctls:
            self.logger.info(
                f"Cleaning up {len(tmp_ctls)} residual tmp*.ctl files from previous runs"
            )
            for tmp_ctl in tmp_ctls:
                try:
                    os.unlink(tmp_ctl)
                    self.logger.debug(f"Removed residual file: {tmp_ctl}")
                except Exception as e:
                    self.logger.warning(
                        f"Failed to remove residual file {tmp_ctl}: {e}"
                    )

    def _manual_protein_fix(self) -> None:
        """
        手动修复蛋白质序列的 Hessian 计算 (仅 PAML 4.10.8)

        PAML 4.10.8 处理蛋白质序列时，mcmctree 生成的 tmp*.ctl 文件
        经常会丢失 aaRatefile 参数，导致 codeml 运行失败。
        此工作流扫描并修复所有 tmp*.ctl 文件，然后手动运行 codeml。

        注意: PAML 4.8 不存在此问题，不需要手动修复。
        """
        import glob

        self.logger.info("Starting manual protein fix workflow for Hessian calculation")

        # 查找所有 tmp*.ctl 文件
        tmp_ctls = glob.glob(str(self.work_dir / "tmp*.ctl"))

        if not tmp_ctls:
            self.logger.warning("No tmp*.ctl files found for protein fix")
            return

        self.logger.info(f"Found {len(tmp_ctls)} tmp*.ctl files to process")

        fixed_count = 0
        for ctl in tmp_ctls:
            ctl_path = Path(ctl)

            # 读取控制文件
            with open(ctl_path, "r") as f:
                content = f.read()

            # 检查是否缺失 aaRatefile
            if "aaRatefile" not in content and "model = 2" in content:
                # 添加 aaRatefile（使用 WAG 替换矩阵作为默认）
                content = content.replace(
                    "model = 2", "model = 2\naaRatefile = wag.dat"
                )

                with safe_writer(ctl_path) as f:
                    f.write(content)

                fixed_count += 1
                self.logger.debug(f"Fixed aaRatefile in {ctl_path.name}")

            # 运行 codeml
            cmd = [self._get_executable("codeml"), ctl_path.name]
            runner = ProcessRunner(
                cwd=self.work_dir,
                timeout=self.common_config.timeout or 3600,
                env=self._get_paml_env(),
            )
            result = runner.run(cmd)

            if result.returncode != 0:
                self.logger.warning(
                    f"codeml failed for {ctl_path.name}: {result.stderr}"
                )

        self.logger.info(f"Fixed {fixed_count} control files, executed codeml for all")

        # 收集所有 rst1 文件
        rst1_files = glob.glob(str(self.work_dir / "rst1*"))

        if rst1_files:
            # 合并 rst1 文件为 in.BV
            with safe_writer(self.work_dir / "in.BV") as out:
                for rst1 in sorted(rst1_files):
                    with open(rst1, "r") as f:
                        out.write(f.read())

            # 复制为 out.BV（供后续阶段使用）
            shutil.copy(self.work_dir / "in.BV", self.work_dir / "out.BV")

            self.logger.info(f"Merged {len(rst1_files)} rst1 files into in.BV/out.BV")
        else:
            self.logger.error("No rst1 files generated after protein fix")
            raise ExecutionError(
                "Protein Hessian calculation failed: no rst1 files generated. "
                "This may indicate a problem with the protein alignment or PAML installation."
            )

    def _run_mcmc_stage(self) -> None:
        """运行 MCMC 采样阶段 (usedata = 2)

        PAML 4.10.8 approximate likelihood MCMC workflow:
        - Uses pre-computed in.BV file (branch lengths, gradient, Hessian)
        - Runs MCMC sampling for divergence time estimation
        - Generates FigTree.tre output file with dated tree
        """
        for run_id in range(1, self.config.num_runs + 1):
            run_dir = self.work_dir / f"run{run_id}"
            run_dir.mkdir(exist_ok=True)

            shutil.copy(self._input_files["tree"], run_dir / "mcmctree.tree")
            shutil.copy(self._input_files["seq"], run_dir / "mcmctree.phy")
            shutil.copy(self.work_dir / "in.BV", run_dir / "in.BV")

            ctl_file = run_dir / "mcmctree.ctl"
            self._generate_control_file(ctl_file, usedata=2, run_id=run_id)

            version_info = "unknown"
            mcmctree_exe = self._get_executable("mcmctree")
            runner = ProcessRunner(
                cwd=run_dir,
                timeout=self.common_config.timeout or 7200,
                env=self._get_paml_env(),
            )
            try:
                version_result = runner.get_version(mcmctree_exe, "--version")
                if version_result:
                    version_info = version_result
            except Exception as e:
                # 版本探测失败只影响留档信息，但原因仍需可查（C-8）
                self.logger.debug(f"Could not probe mcmctree version: {e}")

            figtree = run_dir / "FigTree.tre"
            self.logger.execution_provenance(
                tool_name="mcmctree",
                command=f"{mcmctree_exe} {ctl_file} (usedata=2)",
                version=version_info,
                environment={
                    "PAML_version": version_info,
                    "clock_model": str(self.config.clock),
                    "num_runs": str(self.config.num_runs),
                    "work_dir": str(run_dir),
                },
                artifacts={
                    "input_tree": str(self._input_files["tree"]),
                    "input_seq": str(self._input_files["seq"]),
                    "control_file": str(ctl_file),
                    "inBV": str(run_dir / "in.BV"),
                    "output_figtree": str(figtree) if figtree.exists() else None,
                },
            )

            self.logger.info(
                f"Starting MCMC run {run_id}/{self.config.num_runs} (timeout: 7200s)"
            )

            result = runner.run([mcmctree_exe, ctl_file.name])

            if result.returncode != 0:
                raise ExecutionError(
                    f"MCMC run {run_id} failed (returncode={result.returncode}): {result.stderr}"
                )

            figtree = run_dir / "FigTree.tre"
            if figtree.exists():
                self.logger.info(
                    f"MCMC run {run_id} completed, generated FigTree.tre ({figtree.stat().st_size} bytes)"
                )
            else:
                self.logger.warning(f"MCMC run {run_id} did not generate FigTree.tre")

        self.logger.success(f"Completed {self.config.num_runs} MCMC run(s)")

    def _generate_control_file(
        self, ctl_path: Path, usedata: int, run_id: int = 1
    ) -> None:
        """生成 MCMCTree 控制文件

        根据 PAML 版本生成兼容的控制文件:
        - 4.10.8+: 使用简化的 rgene_gamma/sigma2_gamma 格式，移除 finetune
        - 4.8: 保留旧格式以兼容旧版本
        """
        self._cleanup_residual_tmp_ctl_files()

        alpha_mu = self.config.rate_alpha

        is_paml_48 = self._is_paml_48()
        is_paml_4108_or_later = self._is_paml_4108_or_later()

        # rgene_gamma 的 beta: 若阶段 0 完成了速率估算，使用估算速率
        # (Gamma 均值 = rate * rate_scale)；否则使用基于 base_rate 的默认值
        if self._calculated_beta:
            beta_mu = self._calculated_beta
            self.logger.debug(
                f"Using calculated beta from rate estimation: beta_mu={beta_mu:.4f}"
            )
        else:
            time_scale_factor = 1000.0
            base_rate = 0.1

            if hasattr(self.config, "time_unit") and self.config.time_unit == "Ga":
                adjusted_rate = base_rate * time_scale_factor
                beta_mu = alpha_mu / adjusted_rate
                self.logger.debug(
                    f"Adjusted rgene_gamma for Ga time scale: beta_mu={beta_mu:.4f}"
                )
            else:
                beta_mu = alpha_mu / base_rate

        # RootAge 处理
        # MCMCTree 要求在控制文件中设置 RootAge，或者在树文件中设置根校准
        # 但是 MCMCTree 无法解析树文件中的根节点校准标注
        # 所以我们始终在控制文件中设置 RootAge
        # 注意: MCMCTree 要求根节点年龄约束必须是边界约束（bounds），而不是固定年龄
        rootage_line = ""

        # 查找根节点校准
        root_calibration = None
        for cal in self._calibrations:
            if cal.is_root_node:
                root_calibration = cal
                break

        if self._root_annotation_in_tree:
            # PAML 4.10+ 的根约束已直接写入树文件，控制文件中保留 RootAge
            # 会导致 "Only bounds for the root age are implemented" 错误。
            rootage_line = ""
        elif root_calibration:
            # 根据根节点校准的约束类型，设置合适的 RootAge 值。
            # 单位注意：约束年龄（fixed_age/max_age）以 Ma 存储，而 MCMCTree
            # 控制文件中的 RootAge 以 Ga 为单位，写入时必须 ÷1000 转 Ga，
            # 否则 RootAge 会比树内校准大 1000×。
            from ..models.constraints import (
                FixedAgeConstraint,
                MaximumAgeConstraint,
                SoftBoundsConstraint,
                UniformAgeConstraint,
            )

            constraint = root_calibration.age_constraint
            if isinstance(constraint, FixedAgeConstraint):
                # 固定年龄：转换为窄边界约束
                # 使用 '<age' 格式，MCMCTree 会将其解释为上界
                rootage_line = f"RootAge = '<{constraint.fixed_age / 1000.0}'"
            elif isinstance(constraint, UniformAgeConstraint):
                # 均匀分布：使用最大年龄作为上界（Ma→Ga 需 ÷1000）
                rootage_line = f"RootAge = '<{constraint.max_age / 1000.0}'"
            elif isinstance(constraint, MaximumAgeConstraint):
                # 最大年龄：使用 '<max_age>' 格式（Ma→Ga 需 ÷1000）
                rootage_line = f"RootAge = '<{constraint.max_age / 1000.0}'"
            elif isinstance(constraint, SoftBoundsConstraint):
                # 软边界：使用最大年龄作为上界（Ma→Ga 需 ÷1000）
                rootage_line = f"RootAge = '<{constraint.max_age / 1000.0}'"
            else:
                # 其他类型：使用默认值
                rootage_line = f"RootAge = '<{DEFAULT_ROOT_AGE_GA}'"
                self.logger.warning(
                    f"Unsupported root constraint type: {type(constraint).__name__}. "
                    f"Using default RootAge='<{DEFAULT_ROOT_AGE_GA}'."
                )
        else:
            # 没有根节点校准：使用配置中的 root_age 或默认值
            if self.config.root_age:
                rootage_line = f"RootAge = '<{self.config.root_age / 1000.0}'"
            else:
                rootage_line = f"RootAge = '<{DEFAULT_ROOT_AGE_GA}'"
                # 兜底上界的口径说明（审阅报告 C-9）：10 Ga 故意比
                # constraints.MAX_GEOLOGICAL_AGE_MA（4.6 Ga）更宽，目的是让这一硬
                # 上界对所有生物学相关年龄都不构成约束，从而不干扰 time 超先；
                # 代价是根年龄完全由分子钟 + BDparas 先验决定。
                self.logger.warning(
                    f"No root calibration provided. Using default RootAge='<{DEFAULT_ROOT_AGE_GA}' "
                    f"({DEFAULT_ROOT_AGE_GA * 1000:.0f} Ma). This bound is deliberately wider "
                    "than the physical upper limit used by the constraint layer "
                    "(4600 Ma) so that it never truncates the time prior; the price is that "
                    "the ROOT AGE IS THEN DETERMINED ENTIRELY by the molecular clock and the "
                    "BDparas birth-death prior, with no fossil information at all. "
                    "Please provide --root-age (or a root calibration) for meaningful results."
                )

        # 自动模型选择：'auto' 对 DNA 默认使用 GTR+G (7)，蛋白质使用 2
        model_value: Union[str, int] = self.config.model
        if model_value == "auto":
            model_value = 2 if self._seqtype == 2 else 7

        if self._seqtype == 2:
            model_line = f"model = {model_value}\naaRatefile = wag.dat"
        else:
            model_line = f"model = {model_value}"

        # rgene_gamma 和 sigma2_gamma 格式根据版本调整
        # PAML 4.8 和 4.10.8+ 都使用 2 参数格式: shape scale
        # 但 4.8 可能需要第三个参数作为 alpha 默认值
        if is_paml_48:
            # PAML 4.8 使用 2 参数格式: rgene_gamma = alpha beta
            # sigma2_gamma = alpha beta
            rgene_line = f"rgene_gamma  = {alpha_mu:.4f} {beta_mu:.4f}"
            sigma_line = "sigma2_gamma = 1 4.5"
        elif is_paml_4108_or_later:
            # PAML 4.10.8+ 使用 2 参数格式
            rgene_line = f"rgene_gamma  = {alpha_mu:.4f} {beta_mu:.4f}"
            sigma_line = "sigma2_gamma = 1 10"
        else:
            # 其他版本使用 4 参数格式
            rgene_line = f"rgene_gamma  = {alpha_mu:.4f} {beta_mu:.4f} 1.0 0"
            sigma_line = "sigma2_gamma = 1 10 1 0"

        # BDparas: 4.10.8+ 需要指定构造类型 (m 或 c)
        bd_construction = getattr(self.config, "bd_construction", "c")
        if is_paml_4108_or_later:
            bdparas_line = f"BDparas  = 1 1 0.1 {bd_construction}"
        else:
            bdparas_line = "BDparas  = 1 1 0.1"

        # 检查点设置 (4.10.8+)
        checkpoint_line = ""
        if is_paml_4108_or_later and getattr(self.config, "checkpoint_enabled", False):
            checkpoint_mode = getattr(self.config, "checkpoint_mode", 1)
            checkpoint_prob = getattr(self.config, "checkpoint_prob", 0.01)
            checkpoint_file = getattr(self.config, "checkpoint_file", "mcmctree.ckpt1")
            checkpoint_line = (
                f"checkpoint = {checkpoint_mode} {checkpoint_prob} {checkpoint_file}"
            )

        # 复制基因设置
        duplication_line = ""
        if is_paml_4108_or_later and getattr(self.config, "duplication", 0) == 1:
            duplication_line = "duplication = 1"

        # finetune 在 4.10.8+ 已弃用
        finetune_line = ""
        if is_paml_48:
            finetune_line = "finetune = 1: 0.1 0.1 0.1 0.1 0.1 0.1"

        # usedata 格式: 4.10.8+ 使用 "usedata = 2 in.BV" 格式
        # 注意: usedata = 2 和 in.BV 之间有空格, 不是 usedata = 2  in.BV
        usedata_str = f"usedata = {usedata}"
        if usedata == 2:
            usedata_str += " in.BV"

        # 种子策略：优先使用全局配置的 seed（保证可重现）；
        # PAML 要求随机数种子为正奇数，多次独立运行通过 run_id 偏移保证链间种子不同。
        # 未配置时退回负种子（PAML 负种子表示由系统时钟随机生成）。
        if self.common_config.seed is not None:
            mcmc_seed = (int(self.common_config.seed) + run_id - 1) | 1
        else:
            mcmc_seed = -(run_id * 12345)

        ctl_content = f"""seed = {mcmc_seed}
seqfile  = mcmctree.phy
treefile = mcmctree.tree
outfile  = out.txt
mcmcfile = mcmc.txt

ndata    = {self.config.ndata}
seqtype  = {self._seqtype}
{usedata_str}
clock    = {self.config.clock}

{model_line}

alpha    = 0.5
ncatG    = 4
cleandata = {self.config.clean_data}

{rootage_line}

{bdparas_line}
kappa_gamma  = 6 2
alpha_gamma  = 1 1

{rgene_line}
{sigma_line}

{finetune_line}

{checkpoint_line}

{duplication_line}

print    = 1
burnin   = {self.config.burnin}
sampfreq = {self.config.sampfreq}
nsample  = {self.config.nsample}
"""

        with safe_writer(ctl_path) as f:
            f.write(ctl_content)

    def _is_paml_48(self) -> bool:
        """检测是否为 PAML 4.8 版本"""
        version = getattr(self.config, "paml_version", "4.10.8")
        try:
            parts = version.split(".")
            major = int(parts[0]) if len(parts) > 0 else 4
            minor_str = parts[1] if len(parts) > 1 else "10"
            minor = int(minor_str.split("a")[0])
            return major == 4 and minor == 8
        except (ValueError, IndexError):
            return False

    def _is_paml_4108_or_later(self) -> bool:
        """检测是否为 PAML 4.10.8 或更高版本"""
        version = getattr(self.config, "paml_version", "4.10.8")
        try:
            parts = version.split(".")
            major = int(parts[0]) if len(parts) > 0 else 4
            minor_str = parts[1] if len(parts) > 1 else "10"
            minor = int(minor_str.split("a")[0])
            patch_str = parts[2] if len(parts) > 2 else "8"
            patch = int(patch_str.split("a")[0])

            if major < 4:
                return False
            if major > 4:
                return True
            if major == 4 and minor < 10:
                return False
            if major == 4 and minor > 10:
                return True
            if major == 4 and minor == 10 and patch >= 8:
                return True
            return False
        except (ValueError, IndexError):
            return False

    def _should_skip_rate_estimation(self) -> bool:
        """检测是否应该跳过速率估算阶段（阶段 0）

        由配置项 ``skip_rate_estimation`` 控制（默认 False=执行估算）。
        PAML 4.8 与 4.10.8 的 baseml/codeml 均支持 clock=1 速率估算并输出
        "Substitution rate is per time unit"（4.9i 曾因 bug 缺失该输出，
        4.9j 已修复，详见 PAML changelog 与使用文档）。
        """
        return bool(getattr(self.config, "skip_rate_estimation", False))

    def parse_results(self) -> DatingResult:
        """解析 MCMCTree 输出"""
        # 收集所有运行的结果
        all_ages: Dict[str, List[Dict[str, float]]] = {}
        runs_with_figtree = 0

        for run_id in range(1, self.config.num_runs + 1):
            run_dir = self.work_dir / f"run{run_id}"
            figtree_file = run_dir / "FigTree.tre"

            if figtree_file.exists():
                runs_with_figtree += 1
                ages = self._parse_figtree(figtree_file)
                for name, age in ages.items():
                    if name not in all_ages:
                        all_ages[name] = []
                    all_ages[name].append(age)

        # 计算汇总统计（多链合并口径与 ci_type 严格对应，见 _summarize_chain_ages）
        node_ages = {}
        parsed_from_annotation = False
        interval_kinds: Dict[str, bool] = {}
        divergent_nodes: Dict[str, List[Tuple[int, int]]] = {}
        ga_to_ma = 1000  # FigTree.tre 以 Ga 为单位，PhyloDater 内部统一用 Ma

        for name, age_list in all_ages.items():
            summary = _summarize_chain_ages(age_list)
            if summary is None:
                continue
            parsed_from_annotation = True

            if summary["n_chains"] > 1 and summary["disjoint_chain_pairs"]:
                divergent_nodes[name] = summary["disjoint_chain_pairs"]
            if len(age_list) < self.config.num_runs:
                self.logger.warning(
                    f"MCMCTree node '{name}': only {len(age_list)} of "
                    f"{self.config.num_runs} configured runs reported an age; the "
                    "summary below is based on those runs only."
                )
            interval_kinds[name] = summary["is_envelope"]

            node_ages[name] = NodeAgeEstimate(
                mean_age=summary["mean"] * ga_to_ma,
                median_age=summary["median"] * ga_to_ma,
                ci_lower=summary["ci_lower"] * ga_to_ma,
                ci_upper=summary["ci_upper"] * ga_to_ma,
                ci_type=summary["ci_type"],
            )

        if divergent_nodes:
            detail = "; ".join(
                f"{node} (chains {pairs})"
                for node, pairs in sorted(divergent_nodes.items())
            )
            self.logger.error(
                "MCMCTree multi-chain conflict: the 95% intervals reported by "
                f"different chains do not overlap at {len(divergent_nodes)} node(s): "
                f"{detail}. Non-overlapping chain posteriors are the classic sign of "
                "non-convergence; the pooled interval below cannot be interpreted as "
                "a posterior summary. Increase nsample/burnin or check the priors."
            )
        if any(interval_kinds.values()):
            self.logger.warning(
                "MCMCTree: for nodes sampled by >=2 chains, ci_lower/ci_upper are the "
                "OUTER ENVELOPE of each chain's 95% interval (min of lowers, max of "
                "uppers), NOT a 95% HPD — coverage of an envelope is not defined and "
                "grows with the number of chains. Those nodes are therefore labelled "
                "ci_type=RANGE. Single-chain nodes keep MCMCTree's own HPD (HPD95). "
                "Pooling per-sample node draws would be required to report a genuine "
                "multi-chain HPD; FigTree.tre only stores per-chain summaries."
            )

        if not node_ages and runs_with_figtree:
            self.logger.debug(
                f"MCMCTree: {runs_with_figtree} FigTree.tre file(s) found but no node "
                "age could be extracted from them (chronogram fallback attempted below)."
            )

        # ------------------------------------------------------------------
        # chronogram 回退：若基于 height= 注解 / NEXUS 校准标注的解析未能提取
        # 到任何年龄（例如 FigTree.tre 被截断、输出为 NEXUS 而非 Newick、或
        # 节点名与校准名不匹配），则对时间树做全树 chronogram 解码：以
        # root-to-node 距离（Ga→Ma，×1000）作为节点年龄，避免分析成功却被
        # 报告为 "Nodes dated: 0"。
        #
        # mcmctree 输出树的叶节点为短名（prepare_inputs 用 name_mapping 将全名
        # 压缩为 ACCESSION 短名），而校准点 mrca_leaf_pair 使用全名。定位 MRCA
        # 时先用 name_mapping 把全名转成短名再在树中查找；若映射不可用则直接
        # 用全名尝试。节点键优先取校准点名，其次节点 label，最后 node_<id>。
        # ------------------------------------------------------------------
        if not parsed_from_annotation:
            fallback_tree = None
            try:
                first_run_figtree = self.work_dir / "run1" / "FigTree.tre"
                if first_run_figtree.exists():
                    from dendropy import Tree as _DendropyTree

                    figtree_txt = first_run_figtree.read_text().strip()
                    # 兼容 PAML 4.8 NEXUS 与 4.10.8+ 标准 Newick 两种格式
                    tree_data = figtree_txt
                    nexus_hit = re.search(
                        r"UTREE\s+\d+\s*=\s*(.+);", figtree_txt, re.DOTALL
                    )
                    if nexus_hit:
                        tree_data = nexus_hit.group(1).strip()
                    fallback_tree = _DendropyTree.get(
                        data=tree_data, schema="newick", preserve_underscores=True
                    )
            except Exception as e:
                self.logger.debug(
                    f"MCMCTree chronogram fallback could not load FigTree.tre: {e}"
                )

            if fallback_tree is not None:
                fallback_ages = self._decode_mcmctree_chronogram(fallback_tree)
                if fallback_ages:
                    self.logger.info(
                        "MCMCTree: annotation-based parsing extracted no ages; "
                        f"fell back to chronogram distance decoding ({len(fallback_ages)} nodes)"
                    )
                    node_ages.update(fallback_ages)

        if not node_ages:
            self.logger.warning(
                "MCMCTree: dated tree produced but no node ages could be extracted. "
                "FigTree.tre may be truncated or tip names may not match calibrations."
            )

        # 读取定年树（使用第一次运行的结果）
        first_run_dir = self.work_dir / "run1"
        figtree_file = first_run_dir / "FigTree.tre"
        dated_tree_newick = ""
        if figtree_file.exists():
            file_size = figtree_file.stat().st_size
            # 低于 64 字节的 Newick/NEXUS 树几乎不可能包含有效拓扑（单行最简
            # 树也需远超此数），才提示可能截断。更小的阈值会误报合法的小树
            # （例如 8 叶 NEXUS 树约 570 字节）。
            if file_size < 64:
                self.logger.warning(
                    f"Generated FigTree.tre is suspiciously small ({file_size} bytes). "
                    "Tree file may have been truncated or generation failed."
                )
            with open(figtree_file, "r") as f:
                dated_tree_newick = f.read().strip()

            # PAML 4.8 / 部分输出使用非标准 UTREE 关键字，将其转换为标准
            # NEXUS 的 TREE 关键字，以便 DendroPy/Bio.Phylo 等标准库解析。
            if (
                dated_tree_newick.upper().startswith("#NEXUS")
                and "UTREE" in dated_tree_newick.upper()
            ):
                dated_tree_newick = self._normalize_mcmctree_nexus(dated_tree_newick)
                self.logger.info("Normalized MCMCTree NEXUS output: UTREE -> TREE")

            # PAML 输出树使用与 PHYLIP 序列一致的短名（如 GB_GCA_020632575_1），
            # 这里把叶节点短名恢复为原始全名，使输出树与输入树标签一致。
            dated_tree_newick = self._restore_leaf_names(dated_tree_newick)

        is_converged = self._check_convergence()
        metadata: Dict[str, Any] = {
            "num_runs": self.config.num_runs,
            "clock": self.config.clock,
            "seqtype": self._seqtype,
            "runs_with_figtree": runs_with_figtree,
        }
        if self._convergence_diagnostics:
            # PSRF / ESS / 实际生效阈值随结果一并留档（审阅报告 B-4：阈值不能是
            # 只藏在代码里的魔法数字，报告必须能复述判定依据）
            metadata["convergence_diagnostics"] = dict(self._convergence_diagnostics)
        if interval_kinds:
            metadata["node_age_interval_kinds"] = {
                name: (
                    _ENVELOPE_INTERVAL_DESCRIPTION
                    if is_envelope
                    else _SINGLE_CHAIN_INTERVAL_DESCRIPTION
                )
                for name, is_envelope in interval_kinds.items()
            }
        if divergent_nodes:
            metadata["chain_interval_conflicts"] = {
                name: [list(pair) for pair in pairs]
                for name, pairs in sorted(divergent_nodes.items())
            }

        result = DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=0.0,
            is_converged=is_converged,
            metadata=metadata,
        )

        self.logger.success(f"Parsed MCMCTree results: {len(node_ages)} node ages")
        return result

    def _decode_mcmctree_chronogram(
        self, chronogram: "dendropy.Tree"
    ) -> Dict[str, "NodeAgeEstimate"]:
        """对 mcmctree 输出的 chronogram 做全树解码。

        mcmctree 输出树为 chronogram（年龄编码在分支长度中，单位为 Ga，根=0
        向叶递增）且叶节点使用短名（``prepare_inputs`` 通过
        ``self._name_mapping`` 将全名压缩为 ACCESSION 短名）。本方法遍历所有
        节点，用 ``distance_from_root() × 1000`` 得到节点年龄（Ga→Ma），并优
        先以校准点名作键：先用 name_mapping 把校准 mrca_leaf_pair 全名转为短名
        后在树中定位 MRCA；若某节点与某校准 MRCA 重合则取该校准点名。
        未匹配到校准名的节点依次回退到 ``node.label`` 与 ``node_<id>``。
        """
        from ..models.results import CIType, NodeAgeEstimate

        node_ages: Dict[str, "NodeAgeEstimate"] = {}
        if chronogram is None:
            return node_ages

        # 全名 → 短名 映射（与 prepare_inputs 喂给 mcmctree 的序列名一致）。
        name_mapping = getattr(self, "_name_mapping", None) or {}

        # mcmctree 输出树对 dendropy 表现为无根（is_rooted 未置位），
        # 影响 MRCA 计算与 distance_from_root 语义。在 taxon_namespace 层面
        # 标注为有根，与 mcmctree 的实际有根输出一致。
        try:
            if chronogram.is_rooted is None or chronogram.is_rooted is False:
                chronogram.is_rooted = True
        except Exception as e:
            # C-8：吞掉的异常等于把失败伪装成「没有这回事」，至少留一条可查记录
            self.logger.debug(
                f"MCMCTree chronogram fallback: could not mark tree as rooted: {e}"
            )

        # node id → 校准点名：在树中定位各校准的 MRCA 节点。
        mrca_nodeid_to_calname: Dict[int, str] = {}
        for cal in self._calibrations:
            if not cal.name or not cal.mrca_leaf_pair:
                continue
            taxa = []
            for full in cal.mrca_leaf_pair:
                short = name_mapping.get(full, full)
                taxon = chronogram.taxon_namespace.get_taxon(short)
                if taxon is None:
                    # 短名找不到时回退尝试全名（兼容未改名的情形）
                    taxon = chronogram.taxon_namespace.get_taxon(full)
                if taxon is None:
                    break
                taxa.append(taxon)
            if len(taxa) == len(cal.mrca_leaf_pair):
                try:
                    mrca = chronogram.mrca(taxa=taxa)
                    if mrca is not None:
                        mrca_nodeid_to_calname[id(mrca)] = cal.name
                except Exception as e:
                    # 定位失败 → 该节点不会以校准点名入表，而是退化为 node_<id>；
                    # 记录原因，避免「报告里看不到校准名」时无从追查（C-8）
                    self.logger.debug(
                        f"MCMCTree chronogram fallback: cannot locate MRCA for "
                        f"calibration '{cal.name}': {e}"
                    )

        ga_to_ma = 1000.0
        for node in chronogram.preorder_node_iter():
            try:
                age = float(node.distance_from_root()) * ga_to_ma
            except Exception:
                continue

            # 优先：该节点是某校准的 MRCA → 用校准点名（内部节点，唯一）
            cal_name = mrca_nodeid_to_calname.get(id(node))
            if cal_name and cal_name not in node_ages:
                node_ages[cal_name] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)
                continue

            # 叶节点：以其短名 taxon 标签作键（与 MRCA 点名互不冲突）。
            # 不把校准点名赋给叶节点，否则会与 1085 行的 MRCA 节点 key 冲突、
            # 导致叶节点被去重丢弃。
            tip_taxon = node.taxon
            node_label = node.label or ""
            if tip_taxon is not None and not node.child_nodes():
                key = tip_taxon.label or node_label
            else:
                # 内部分支：优先节点 label，再 node_<id>
                key = node_label if node_label else f"node_{id(node)}"

            if key and key not in node_ages:
                node_ages[key] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)

        return node_ages

    def _parse_figtree(self, figtree_path: Path) -> Dict[str, Dict[str, float]]:
        """解析 FigTree.tre 文件，通过 MRCA 拓扑匹配校准点

        同时支持：
          - PAML 4.10.8+ 的标准 Newick（节点注解为 ``[&95%HPD={lower,upper}]``）
          - PAML 4.8 的 NEXUS/UTREE 输出（注解为 ``[&height=...,height_95%_HPD={...}]``）

        当注解只有 HPD 区间时，取其中点作为 mean；解析结果单位为 Ga，外层会
        乘以 1000 转回 Ma。
        """
        ages: Dict[str, Dict[str, float]] = {}

        with open(figtree_path, "r") as f:
            content = f.read().strip()

        # ------------------------------------------------------------------
        # 路径 A: PAML 4.8 NEXUS 格式 (#NEXUS  BEGIN TREES;  UTREE 1 = ...)
        # ------------------------------------------------------------------
        if content.startswith("#NEXUS") or "UTREE" in content.split("\n", 3)[:3]:
            name_mapping = getattr(self, "_name_mapping", None) or {}
            return _parse_figtree_nexus(
                content,
                self._calibrations,
                logger=self.logger,
                name_mapping=name_mapping,
            )

        # ------------------------------------------------------------------
        # 路径 B: PAML 4.10.8+ 标准 Newick
        # ------------------------------------------------------------------
        try:
            from dendropy import Tree as _DendropyTree

            tree = _DendropyTree.get(
                data=content, schema="newick", preserve_underscores=True
            )
            _ensure_rooted(tree)
        except Exception as e:
            self.logger.warning(f"Could not parse FigTree.tre with dendropy: {e}")
            return ages

        name_mapping = getattr(self, "_name_mapping", None) or {}
        for cal in self._calibrations:
            if not cal.name:
                continue
            try:
                node = self._locate_calibration_node(
                    tree, cal, name_mapping=name_mapping
                )
                if node is None:
                    continue
                age = _annotation_to_age(node)
                if age:
                    ages[cal.name] = age
            except Exception as e:
                self.logger.debug(
                    f"MCMCTree: could not extract age for {cal.name}: {e}"
                )

        # 补充未命名的内部节点（非校准节点），使 comparison_table 与其他方法一致
        try:
            calibrated_node_ids = set()
            for cal in self._calibrations:
                if not cal.name:
                    continue
                try:
                    node = self._locate_calibration_node(
                        tree, cal, name_mapping=name_mapping
                    )
                    if node is not None:
                        calibrated_node_ids.add(id(node))
                except Exception as e:
                    # 定位失败只影响「哪些内部节点算已校准」的去重，但留下痕迹
                    self.logger.debug(
                        f"MCMCTree: could not locate calibrated node for {cal.name}: {e}"
                    )

            internal_idx = 1
            for node in tree.preorder_node_iter():
                if node.is_leaf():
                    continue
                if id(node) in calibrated_node_ids:
                    continue
                age = _annotation_to_age(node)
                if age:
                    while f"internal_node_{internal_idx}" in ages:
                        internal_idx += 1
                    ages[f"internal_node_{internal_idx}"] = age
                    internal_idx += 1
        except Exception as e:
            self.logger.debug(
                f"MCMCTree: could not extract unnamed internal node ages: {e}"
            )

        return ages

    def _locate_calibration_node(
        self,
        tree: "dendropy.Tree",
        cal: CalibrationPoint,
        name_mapping: Optional[Dict[str, str]] = None,
    ) -> Optional["dendropy.Node"]:
        """在校准树中定位校准点对应的节点"""
        if name_mapping is None:
            name_mapping = {}
        if cal.is_root_node:
            return tree.seed_node
        if not cal.mrca_leaf_pair or len(cal.mrca_leaf_pair) < 2:
            return None
        taxa = []
        for tip in cal.mrca_leaf_pair:
            short = name_mapping.get(tip, tip)
            taxon = tree.taxon_namespace.get_taxon(short)
            if taxon is None:
                break
            taxa.append(taxon)
        if len(taxa) != len(cal.mrca_leaf_pair):
            return None
        return tree.mrca(taxa=taxa)

    def _check_convergence(self) -> Optional[bool]:
        """检查 MCMC 收敛性：Gelman–Rubin PSRF + Geyer ESS，三态返回（A-1/A-2/B-4）。

        返回语义（B-4 要求把「布尔判定」与「诊断量」分开）：
          True  — 每条随机参数列都满足 PSRF <= psrf_threshold 且 ESS >= ess_threshold
          False — 未达阈值、链含 NaN、存在真正停滞（零方差）的随机参数列，
                  或 mcmc.txt 缺失/无法解析（无法核验即按未收敛处理，fail-closed）
          None  — 无法诊断：num_runs < 2、可用链不足 2 条、各链列不可对齐、numpy 不可用

        与旧实现的差异：
        1. A-1：MCMCTree 的 mcmc.txt 第一行是**非数值表头**（由 ``collectx()`` 写出，
           见 mcmctree.c 4.8:1724-1741 / 4.10.8:2045-2110）。旧代码
           ``np.loadtxt(mcmc_file)`` 在 'Gen' 上抛 ValueError → 第一条链即
           ``return False``，PSRF/ESS 一次都没有被计算过，诊断系统形同装饰。
        2. A-2：迭代号列（``Gen``）不是后验抽样，各链逐元素相同（1000, 2000, …）。
           留在诊断里有两个危害：链间方差恒为 0 使该列 PSRF ≈ 1（等于用一条
           「永远达标」的假参数稀释判定），而它的组内 ESS 极低（单调趋势的自相关
           时间极长），于是旧实现每次正常运行都只会打出 "PSRF OK but ESS may be low"
           然后照样 return True（把诊断失效伪装成「有点小问题」）；阈值收紧后若不
           剔除，任何正常运行都会被这一列判成未收敛。现按表头列名剔除，并对无表头
           的历史/异常输出用「严格单调递增整数」这一形态判据兜底剔除；真正停滞的
           （零方差）随机参数列仍然明确报错。
        3. B-4：旧实现 PSRF=1.9（<2.0）或 ESS 远低于阈值时仍 ``return True``。
           现只在两项阈值同时满足时返回 True；阈值可由配置项
           ``psrf_threshold`` / ``ess_threshold`` 覆盖，实测值与生效阈值一起写入
           ``self._convergence_diagnostics``，供报告打印实际判定依据。
        4. 旧的 ``if len(all_chains) < 2: return True`` 既不可达、方向也是错的
           （链数不足应是「无法诊断」，不是「已收敛」）。

        注（审阅报告 E-5，只记录不改动）：burn-in 样本不会被写进 mcmc.txt
        （上游仅在 ``ir >= 0`` 时打印样本，见 4.10.8 mcmctree.c:4246-4250），
        因此这里无需、也不应再做 burn-in 截断。

        参考文献:
        - Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
          using multiple sequences. Statistical Science, 7(4), 457-472.
        """
        self._convergence_diagnostics = None

        if self.config.num_runs < 2:
            self.logger.info(
                "MCMC convergence diagnostic requires num_runs >= 2; "
                "returning None for single run."
            )
            return None

        try:
            import numpy as np  # noqa: F401  (与模块级导入一致，仅为保留降级路径)
        except ImportError:
            # 无法诊断收敛：返回 None（未知）而非 True（已收敛），
            # 避免在无 numpy 环境下静默误判为「已收敛」。下游对
            # is_converged 为 None 的情况按「未知」处理，不会误判失败。
            self.logger.warning(
                "numpy 不可用，无法执行 MCMC 收敛诊断（ESS/PSRF）。"
                "is_converged 将标记为未知（None），请根据其他指标谨慎评估结果。"
            )
            return None

        psrf_threshold = _positive_config(
            self.config, "psrf_threshold", DEFAULT_PSRF_THRESHOLD
        )
        ess_threshold = _positive_config(
            self.config, "ess_threshold", DEFAULT_ESS_THRESHOLD
        )

        all_chains: List["np.ndarray"] = []
        column_names: Optional[List[str]] = None

        for run_id in range(1, self.config.num_runs + 1):
            run_dir = self.work_dir / f"run{run_id}"
            mcmc_file = run_dir / "mcmc.txt"

            if not mcmc_file.exists():
                self.logger.warning(
                    f"mcmc.txt not found for run {run_id} (expected {mcmc_file}). "
                    "Convergence cannot be verified; treating the run as NOT converged."
                )
                return False

            try:
                samples, names = _read_mcmc_chain(mcmc_file)
            except Exception as e:
                self.logger.error(
                    f"Failed to load mcmc.txt for run {run_id} ({mcmc_file}): {e}. "
                    "MCMC convergence cannot be verified; treating the run as NOT converged."
                )
                return False

            if column_names is None:
                column_names = names
            elif names != column_names:
                self.logger.error(
                    f"mcmc.txt column layout differs between run 1 and run {run_id}: "
                    f"{column_names} vs {names}. PSRF/ESS need aligned parameter "
                    "columns; convergence cannot be diagnosed."
                )
                return None

            all_chains.append(samples)

        if len(all_chains) < 2:
            # 与旧实现相反：链数不足是「无法诊断」，不是「已收敛」
            self.logger.warning(
                f"Only {len(all_chains)} usable MCMC chain(s); convergence cannot be "
                "assessed (is_converged = None)."
            )
            return None

        lengths = sorted({int(chain.shape[0]) for chain in all_chains})
        if len(lengths) > 1:
            self.logger.error(
                f"MCMC chains have different sample counts {lengths}; the Gelman-Rubin "
                "PSRF implemented here assumes equal-length chains. Convergence "
                "cannot be diagnosed (is_converged = None)."
            )
            return None

        # 检查数据完整性和有效性
        for run_id, chain in enumerate(all_chains, 1):
            if np.any(np.isnan(chain)):
                self.logger.error(
                    f"MCMC chain {run_id} contains NaN values. "
                    "This indicates numerical instability or insufficient burn-in."
                )
                self._convergence_diagnostics = {
                    "error": f"chain {run_id} contains NaN values",
                    "n_chains": len(all_chains),
                    "verdict": "not_converged",
                }
                return False

        diagnostics = _diagnose_mcmc_chains(
            all_chains, column_names or [], psrf_threshold, ess_threshold
        )

        self._convergence_diagnostics = {
            "n_chains": len(all_chains),
            "n_samples_per_chain": int(all_chains[0].shape[0]),
            "psrf_threshold": psrf_threshold,
            "ess_threshold": ess_threshold,
            "max_psrf": _finite_or_none(diagnostics["max_psrf"]),
            "min_ess": _finite_or_none(diagnostics["min_ess"]),
            "psrf_by_column": {
                name: _finite_or_none(value)
                for name, value in diagnostics["psrf_by_column"].items()
            },
            "zero_variance_columns": diagnostics["zero_variance_columns"],
            "index_columns_removed": diagnostics["index_columns"],
        }

        # 报告零方差参数（迭代号一类的确定性索引列已在读取/诊断阶段剔除，
        # 因此走到这里的都是真正的随机参数列停滞 → 明确判未收敛）
        if diagnostics["zero_variance_columns"]:
            stuck = diagnostics["zero_variance_columns"]
            self.logger.error(
                f"MCMC convergence check failed: {len(stuck)} stochastic parameter(s) "
                f"have zero within-chain variance in every chain "
                f"(columns: {stuck[:10]}{'...' if len(stuck) > 10 else ''}). "
                "This indicates stuck chains - the MCMC sampler is not exploring the "
                "parameter space properly. Possible causes: burn-in too short, chain "
                "starting too close to boundary, or numerical overflow."
            )
            self._convergence_diagnostics["verdict"] = "not_converged"
            self._convergence_diagnostics["reasons"] = [
                "zero-variance stochastic parameter column(s): " + ", ".join(stuck)
            ]
            return False

        max_psrf = diagnostics["max_psrf"]
        min_ess = diagnostics["min_ess"]
        if max_psrf is None or min_ess is None:
            self.logger.warning(
                "MCMC chains contain no stochastic parameter columns to diagnose "
                f"(index-like columns removed: {diagnostics['index_columns']}); "
                "convergence cannot be assessed (is_converged = None)."
            )
            self._convergence_diagnostics["verdict"] = "not_diagnosable"
            return None

        self.logger.info(
            f"MCMC diagnostics over {len(all_chains)} chains × "
            f"{int(all_chains[0].shape[0])} samples: max PSRF = {max_psrf:.4f} "
            f"(threshold <= {psrf_threshold}), min ESS = {min_ess:.1f} "
            f"(threshold >= {ess_threshold})"
        )

        if max_psrf <= psrf_threshold and min_ess >= ess_threshold:
            self.logger.success(
                f"MCMC chains have converged (max PSRF {max_psrf:.4f} <= "
                f"{psrf_threshold}, min ESS {min_ess:.1f} >= {ess_threshold})"
            )
            self._convergence_diagnostics["verdict"] = "converged"
            self._convergence_diagnostics["reasons"] = []
            return True

        reasons: List[str] = []
        if max_psrf > psrf_threshold:
            reasons.append(
                f"max PSRF {max_psrf:.4f} exceeds threshold {psrf_threshold} "
                "(chains have not mixed to the same distribution)"
            )
        if min_ess < ess_threshold:
            reasons.append(
                f"min ESS {min_ess:.1f} is below threshold {ess_threshold:.0f} "
                "(posterior summaries are not reliable)"
            )
        # 注意日志文案与返回值必须一致：旧实现写「may not have fully converged」
        # 却 return True，下游 is_converged 因此被记为真。
        self.logger.error(
            "MCMC chains have NOT converged: "
            + "; ".join(reasons)
            + f". Diagnostics: {self._convergence_diagnostics}"
        )
        self._convergence_diagnostics["verdict"] = "not_converged"
        self._convergence_diagnostics["reasons"] = reasons
        return False

    def _normalize_mcmctree_nexus(self, nexus_text: str) -> str:
        """将 PAML 非标准 NEXUS 输出转为标准 NEXUS。

        PAML 4.8 及部分 4.10.8 输出使用 UTREE 关键字，且 :/，后常带空格，
        标准 NEXUS 解析器（如 DendroPy）无法识别。本方法：
          1) 将 UTREE 替换为 TREE；
          2) 移除 : 和 , 后的多余空格；
          3) 保持 #NEXUS / BEGIN TREES / END; 结构。
        """
        import re as _re

        # 1. 把 UTREE 关键字替换为 TREE
        text = _re.sub(r"\bUTREE\b", "TREE", nexus_text, flags=_re.IGNORECASE)

        # 2. 移除分支长度/逗号后的空格，使其符合标准 Newick 紧凑格式
        text = _re.sub(r":\s+", ":", text)
        text = _re.sub(r",\s+", ",", text)

        return text

    def _restore_leaf_names(self, nexus_text: str) -> str:
        """把 MCMCTree 输出 NEXUS 树中的短名叶节点恢复为原始全名。

        prepare_inputs 阶段为 PHYLIP 安全把全名压缩为短名（如
        ``GB_GCA_020632575_1``），PAML 输出树仍使用这些短名。本方法利用
        ``_name_mapping`` 建立反向映射，通过 DendroPy 读取 NEXUS、修改 taxon
        标签后写回，保留节点注解（HPD 区间等）。
        """
        name_mapping = getattr(self, "_name_mapping", None)
        if not name_mapping:
            return nexus_text

        reverse_mapping = {short: full for full, short in name_mapping.items()}

        try:
            from dendropy import Tree as _DendropyTree

            tree = _DendropyTree.get(
                data=nexus_text,
                schema="nexus",
                preserve_underscores=True,
            )

            # 更新 taxon_namespace 中所有叶节点标签
            for taxon in tree.taxon_namespace:
                full_name = reverse_mapping.get(taxon.label)
                if full_name:
                    taxon.label = full_name

            # 同时更新树节点上的 taxon 引用标签
            for node in tree.leaf_nodes():
                if node.taxon:
                    full_name = reverse_mapping.get(node.taxon.label)
                    if full_name:
                        node.taxon.label = full_name

            # dendropy 无类型标注：as_string() 运行期返回 str。
            rendered: str = tree.as_string(schema="nexus")
            return rendered
        except ImportError:
            self.logger.warning(
                "DendroPy not available, cannot restore MCMCTree leaf names"
            )
            return nexus_text
        except Exception as e:
            self.logger.warning(f"Failed to restore MCMCTree leaf names: {e}")
            return nexus_text


# 注册适配器
DatingMethodRegistry.register("mcmctree", MCMCTreeMethod)


def _parse_figtree_nexus(
    content: str,
    calibrations: List["CalibrationPoint"],
    logger: Optional[Any] = None,
    name_mapping: Optional[Dict[str, str]] = None,
) -> Dict[str, Dict[str, float]]:
    """Parse PAML's NEXUS-formatted FigTree.tre output.

    同时兼容 PAML 4.8 (``height=...``) 与 PAML 4.10.8+ (``95%HPD={...}``) 两种
    节点注解。UTREE 语句会被提取并清理为 dendropy 可解析的标准 Newick。
    """
    import re as _re

    ages: Dict[str, Dict[str, float]] = {}

    utree_match = _re.search(
        r"(?:UTREE|TREE)\s+\d+\s*=\s*(.+?);", content, _re.DOTALL | _re.IGNORECASE
    )
    if not utree_match:
        return ages

    treetext = utree_match.group(1).strip()
    # PAML 输出在 : 和 , 后加了空格，dendropy newick 解析器需要紧凑格式
    treetext = _re.sub(r":\s+", ":", treetext)
    treetext = _re.sub(r",\s+", ",", treetext)
    if not treetext.endswith(";"):
        treetext += ";"

    try:
        from dendropy import Tree as _DendropyTree

        tree = _DendropyTree.get(
            data=treetext, schema="newick", preserve_underscores=True
        )
        _ensure_rooted(tree)
    except Exception as e:
        if logger is not None:
            logger.warning(f"Could not parse NEXUS FigTree.tre with dendropy: {e}")
        return ages

    if name_mapping is None:
        name_mapping = {}

    for cal in calibrations:
        if not cal.name:
            continue
        try:
            node = None
            if cal.is_root_node:
                node = tree.seed_node
            elif cal.mrca_leaf_pair and len(cal.mrca_leaf_pair) >= 2:
                taxa = []
                for tip in cal.mrca_leaf_pair:
                    short = name_mapping.get(tip, tip)
                    taxon = tree.taxon_namespace.get_taxon(short)
                    if taxon is None:
                        break
                    taxa.append(taxon)
                if len(taxa) == len(cal.mrca_leaf_pair):
                    node = tree.mrca(taxa=taxa)
            if node is None:
                continue
            age = _annotation_to_age(node)
            if age:
                ages[cal.name] = age
        except Exception as e:
            if logger is not None:
                logger.debug(
                    f"MCMCTree NEXUS: could not extract age for {cal.name}: {e}"
                )

    # 补充未命名的内部节点（非校准节点）
    try:
        calibrated_node_ids = set()
        for cal in calibrations:
            if not cal.name:
                continue
            node = None
            if cal.is_root_node:
                node = tree.seed_node
            elif cal.mrca_leaf_pair and len(cal.mrca_leaf_pair) >= 2:
                taxa = []
                for tip in cal.mrca_leaf_pair:
                    short = name_mapping.get(tip, tip)
                    taxon = tree.taxon_namespace.get_taxon(short)
                    if taxon is None:
                        break
                    taxa.append(taxon)
                if len(taxa) == len(cal.mrca_leaf_pair):
                    node = tree.mrca(taxa=taxa)
            if node is not None:
                calibrated_node_ids.add(id(node))

        internal_idx = 1
        for node in tree.preorder_node_iter():
            if node.is_leaf():
                continue
            if id(node) in calibrated_node_ids:
                continue
            age = _annotation_to_age(node)
            if age:
                while f"internal_node_{internal_idx}" in ages:
                    internal_idx += 1
                ages[f"internal_node_{internal_idx}"] = age
                internal_idx += 1
    except Exception as e:
        if logger is not None:
            logger.debug(
                f"MCMCTree NEXUS: could not extract unnamed internal node ages: {e}"
            )

    return ages


def _ensure_rooted(tree: "dendropy.Tree") -> None:
    """将 dendropy 树标记为有根，避免 mrca() 在无根模式下错误返回根节点。"""
    try:
        if tree.is_rooted is None or tree.is_rooted is False:
            tree.is_rooted = True
    except Exception:
        pass


def _annotation_to_age(node: "dendropy.Node") -> Optional[Dict[str, float]]:
    """从 dendropy 节点注解中提取年龄（Ga）和 95% HPD。

    支持：
      - ``height=mean`` + 可选 ``height_95%_HPD={lower,upper}``
      - 仅 ``95%HPD={lower,upper}``（取中点作为 mean）
    """
    anns: Dict[str, Any] = {}
    for ann in node.annotations:
        name = str(ann.name) if ann.name is not None else ""
        anns[name] = ann.value

    # 路径 1: PAML 4.8 / 部分 4.10.8 输出
    if "height" in anns:
        try:
            mean = float(anns["height"])
        except Exception:
            mean = None
        if mean is not None:
            lower = upper = mean
            hpd = anns.get("height_95%_HPD")
            if isinstance(hpd, (list, tuple)) and len(hpd) == 2:
                try:
                    lower = float(hpd[0])
                    upper = float(hpd[1])
                except Exception:
                    pass
            return {"mean": mean, "median": mean, "lower": lower, "upper": upper}

    # 路径 2: PAML 4.10.8+ 仅输出 95% HPD
    hpd = anns.get("95%HPD")
    if isinstance(hpd, (list, tuple)) and len(hpd) == 2:
        try:
            lower = float(hpd[0])
            upper = float(hpd[1])
            mean = (lower + upper) / 2.0
            return {"mean": mean, "median": mean, "lower": lower, "upper": upper}
        except Exception:
            return None

    # 路径 3: PAML 4.8 输出 95%={lower,upper}
    hpd_48 = anns.get("95%")
    if isinstance(hpd_48, (list, tuple)) and len(hpd_48) == 2:
        try:
            lower = float(hpd_48[0])
            upper = float(hpd_48[1])
            mean = (lower + upper) / 2.0
            return {"mean": mean, "median": mean, "lower": lower, "upper": upper}
        except Exception:
            return None

    return None


def _positive_config(config: Any, name: str, default: float) -> float:
    """读取配置中的正数阈值项；缺失/非法时回退默认值并告警由调用方负责。"""
    value = getattr(config, name, None)
    if value is None:
        return float(default)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return float(default)
    if numeric <= 0:
        return float(default)
    return numeric


def _finite_or_none(value: Optional[float]) -> Optional[float]:
    """把 inf/nan 转成 None，保证诊断字典可 JSON 序列化。"""
    if value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(numeric):
        return None
    return numeric


def _all_tokens_numeric(tokens: List[str]) -> bool:
    """判断一组 token 是否全为数值（用于区分 mcmc.txt 的表头行与数据行）。"""
    if not tokens:
        return False
    for token in tokens:
        try:
            float(token)
        except (TypeError, ValueError):
            return False
    return True


def _split_iteration_columns(
    header: List[str],
) -> Tuple[List[int], List[str]]:
    """按 **列名**（而非位置）挑出迭代号列，返回 (保留列下标, 被剔除列名)。

    MCMCTree 的首列是采样迭代号 ``Gen``（``collectx()`` 写表头、采样循环写
    ``ir + 1``）。它不是后验抽样，各链逐元素相同，必须排除在 PSRF/ESS 之外
    （审阅报告 A-2）：它的链间方差为 0（PSRF 恒 ≈ 1，形同一条永远达标的假参数），
    而它的组内 ESS 极低（单调趋势），会把「正常运行」拖成「ESS 不足」。
    """
    keep: List[int] = []
    dropped: List[str] = []
    for idx, name in enumerate(header):
        if name.strip().lower() in _MCMC_ITERATION_COLUMN_NAMES:
            dropped.append(name)
        else:
            keep.append(idx)
    return keep, dropped


def _looks_like_iteration_index(values: "np.ndarray") -> bool:
    """判断一列是否是「严格单调递增的整数索引」（确定性列，非随机参数）。

    用于无表头输入（被截断/老输出）的兜底：MCMCTree 的迭代号列在每个采样点上都
    递增，而任何随机参数列在整个采样序列上单调递增的概率可视为 0。
    """
    if values.size < 2:
        return False
    if not np.all(np.isfinite(values)):
        return False
    if not np.allclose(values, np.round(values)):
        return False
    return bool(np.all(np.diff(values) > 0))


def _read_mcmc_chain(
    mcmc_file: Path,
) -> Tuple["np.ndarray", List[str]]:
    """读取一条 MCMCTree ``mcmc.txt``，返回 (剔除迭代号列后的样本矩阵, 保留列名)。

    真实 mcmc.txt 的第一行是 **非数值表头**::

        Gen\\tt_n5\\tt_n6\\tt_n7\\tmu\\tsigma2\\tlnL

    （见 ``mcmctree.c`` 的 ``collectx()``：4.8:1724-1741、4.10.8:2045-2110；
    InfiniteSites 通路 4.8:1337-1342、4.10.8:1536-1541 同构）。
    ``np.loadtxt`` 默认只把 ``#`` 当注释，因此在第 0 行第 1 列的 'Gen' 上抛
    ``ValueError`` —— 这正是 A-1 里让每一次多链运行都被判「未收敛」的直接原因。

    Raises:
        ValueError: 文件为空、只有表头没有样本、或表头列数与数据列数不一致。
    """
    with open(mcmc_file, "r", encoding="utf-8", errors="ignore") as fh:
        first_line = fh.readline()
        has_data_rows = any(line.strip() for line in fh)

    if not first_line.strip():
        raise ValueError("mcmc.txt is empty (no header, no samples)")

    header = first_line.split()
    if _all_tokens_numeric(header):
        # 没有表头行（异常/被截断的输出）：按位置命名列，稍后用形态判据剔除索引列
        data = np.loadtxt(mcmc_file, ndmin=2)
        header = [f"col{idx + 1}" for idx in range(data.shape[1])]
    else:
        if not has_data_rows:
            raise ValueError(
                f"mcmc.txt contains only a header line ({len(header)} columns); "
                "no MCMC samples were written"
            )
        data = np.loadtxt(mcmc_file, skiprows=1, ndmin=2)
        if data.shape[1] != len(header):
            raise ValueError(
                f"mcmc.txt header has {len(header)} columns but sample rows have "
                f"{data.shape[1]} (header: {' '.join(header)})"
            )

    data = np.atleast_2d(data)
    keep, dropped = _split_iteration_columns(header)
    if not keep:
        raise ValueError(
            f"mcmc.txt has no stochastic parameter columns after removing the "
            f"iteration index column(s) {dropped}"
        )
    return np.ascontiguousarray(data[:, keep]), [header[idx] for idx in keep]


def _diagnose_mcmc_chains(
    chains: List["np.ndarray"],
    column_names: List[str],
    psrf_threshold: float = DEFAULT_PSRF_THRESHOLD,
    ess_threshold: float = DEFAULT_ESS_THRESHOLD,
) -> Dict[str, Any]:
    """对逐参数计算 PSRF 与 ESS，并区分「随机参数停滞」与「确定性索引列」。

    仅做诊断量计算，不做收敛与否的布尔判定（判定在 ``_check_convergence``，
    以便日志文案与返回值一致——见审阅报告 B-4）。

    Args:
        chains: m 条等长、已剔除迭代号列的样本矩阵（形状 n × p）。
        column_names: 保留列的名字（长度 p）。
        psrf_threshold / ess_threshold: 仅用于在返回里回显生效阈值。

    Returns:
        ``{"max_psrf", "min_ess", "psrf_by_column", "zero_variance_columns",
        "index_columns", "psrf_threshold", "ess_threshold"}``；
        无任何随机列时 max_psrf/min_ess 为 None。
    """
    n_params = chains[0].shape[1]
    names = list(column_names[:n_params]) + [
        f"col{i + 1}" for i in range(len(column_names), n_params)
    ]

    psrf_values: List[float] = []
    ess_values: List[float] = []
    psrf_by_column: Dict[str, float] = {}
    zero_variance_columns: List[str] = []
    index_columns: List[str] = []

    for param_idx in range(n_params):
        chains_for_param = [chain[:, param_idx] for chain in chains]
        name = names[param_idx]

        # 确定性索引列（迭代号 / 样本序号）：它不是后验抽样。组内方差很大但完全
        # 由趋势贡献，_geyer_ess 会把这种单调列的 ESS 压到远低于阈值——若不剔除，
        # 每一次正常运行都会因这一列被判「未收敛」（审阅报告 A-2 的另一半危害）。
        if _looks_like_iteration_index(chains_for_param[0]):
            index_columns.append(name)
            continue

        psrf, is_zero_var = _gelman_rubin_psrf(chains_for_param)
        if is_zero_var:
            zero_variance_columns.append(name)
            psrf_by_column[name] = float("inf")
            continue

        psrf_values.append(psrf)
        psrf_by_column[name] = psrf

        # 标准 Geyer (1992) initial positive sequence ESS；
        # 对每条链分别计算，取最小值（最保守估计）。
        ess_values_for_param = [_geyer_ess(c) for c in chains_for_param]
        ess_values.append(min(ess_values_for_param) if ess_values_for_param else 0.0)

    return {
        "max_psrf": max(psrf_values) if psrf_values else None,
        "min_ess": min(ess_values) if ess_values else None,
        "psrf_by_column": psrf_by_column,
        "zero_variance_columns": zero_variance_columns,
        "index_columns": index_columns,
        "psrf_threshold": psrf_threshold,
        "ess_threshold": ess_threshold,
    }


def _summarize_chain_ages(
    age_list: List[Dict[str, float]],
) -> Optional[Dict[str, Any]]:
    """把同一节点在多条 MCMC 链上的摘要合并为一条记录（审阅报告 B-3）。

    口径（并与返回的 ``ci_type`` 严格一一对应）：

    - ``mean``：各链均值的算术平均。
    - ``median``：各链中位数的中位数——**与 mean 使用同一批链**。旧实现 mean 取
      全链、median 只取 ``age_list[0]``（还未必是 run1），同一记录里两个统计量
      来自不同样本集合、不可比。FigTree.tre 只有每链摘要量、无法还原池化样本，
      因此这是对跨链摘要的一致合并，不是真正的池化后验中位数（文档已写明）。
    - 单链（1 条链报告该节点）：区间就是 MCMCTree 自己给出的 95% HPD
      → ``ci_type = HPD95``。
    - 多链（>=2 条）：区间是各链 95% 区间的 **外包络**（min(下界)/max(上界)）。
      其覆盖率随链数增大、没有概率含义，因此标为 ``ci_type = RANGE`` 并在
      ``interval_description`` / 日志 / metadata 里如实说明——绝不冒称 HPD95。
    - ``disjoint_chain_pairs``：任意两条链的 95% 区间完全不相交时的链号对。
      这是典型的不收敛信号，必须暴露而不是被外包络抹平。
    """
    if not age_list:
        return None

    n_chains = len(age_list)
    means = [float(a["mean"]) for a in age_list]
    medians = [float(a.get("median", a["mean"])) for a in age_list]
    lowers = [float(a["lower"]) for a in age_list]
    uppers = [float(a["upper"]) for a in age_list]

    disjoint_pairs: List[Tuple[int, int]] = []
    if n_chains > 1:
        for i in range(n_chains):
            for j in range(i + 1, n_chains):
                if lowers[i] > uppers[j] or lowers[j] > uppers[i]:
                    disjoint_pairs.append((i + 1, j + 1))

    return {
        "mean": sum(means) / n_chains,
        "median": float(np.median(np.asarray(medians, dtype=float))),
        "ci_lower": min(lowers),
        "ci_upper": max(uppers),
        "n_chains": n_chains,
        "ci_type": CIType.HPD95 if n_chains == 1 else CIType.RANGE,
        "is_envelope": n_chains > 1,
        "interval_description": (
            _SINGLE_CHAIN_INTERVAL_DESCRIPTION
            if n_chains == 1
            else _ENVELOPE_INTERVAL_DESCRIPTION
        ),
        "disjoint_chain_pairs": disjoint_pairs,
    }


def _gelman_rubin_psrf(chains: List["np.ndarray"]) -> Tuple[float, bool]:
    """计算单参数的标准 Gelman-Rubin PSRF（Potential Scale Reduction Factor）。

    采用 Gelman et al. (2013, BDA3) 的标准公式：
        W      = mean(s_i^2)                                  # 组内方差
        B/n    = Var(mean_i)  (ddof=1)                         # 组间方差 / n
        var+   = ((n-1)/n) * W + B/n
        PSRF   = sqrt(var+ / W)
    其中 n 为每条链样本数，s_i^2 为第 i 条链方差，mean_i 为第 i 条链均值。

    Args:
        chains: m 条（m>=2）同长度的一维链。

    Returns:
        (psrf, is_zero_variance)。当所有链方差趋零（停滞链）时
        is_zero_variance 为 True，psrf 记为 inf。
    """
    n = chains[0].shape[0]
    variances = [float(np.var(c, ddof=1)) for c in chains]
    if all(v < 1e-15 for v in variances):
        return float("inf"), True
    w = float(np.mean(variances))
    chain_means = [float(np.mean(c)) for c in chains]
    b_over_n = float(np.var(chain_means, ddof=1))
    sigma_hat_sq = ((n - 1) / n) * w + b_over_n
    if w > 1e-15:
        psrf = float(np.sqrt(sigma_hat_sq / w))
    else:
        psrf = 1.0
    return psrf, False


def _geyer_ess(chain: "np.ndarray") -> float:
    """Geyer (1992) initial positive sequence 估计有效样本量（ESS）。

    计算归一化自相关序列 rho_k（rho_0 = 1），构造配对序列
        pair_t = rho_{2t} + rho_{2t+1},
    取其非负且单调不增的前缀求和得到自相关时间
        tau = -1 + 2 * sum_t pair_t,
    则 ESS = n / tau。

    该估计量要求配对序列保持初始正性（initial positive）与单调性
    （initial monotone），因此在首次出现 pair <= 0 或 pair 较前一项增大
    时即停止累加。这与旧实现（从 t=1 起累加、漏掉 rho_1、且无单调/正性
    截断）不同，能给出标准的有效样本量。

    Args:
        chain: 一维 MCMC 样本序列。

    Returns:
        有效样本量（>= 0）。
    """
    n = int(chain.shape[0])
    c_centered = chain - np.mean(chain)
    c_var = float(np.var(c_centered, ddof=1))
    if c_var < 1e-15:
        return 0.0

    # 自相关函数（优先用 FFT 卷积加速，避免长链上的 O(n^2) 纯 Python 开销）。
    try:
        from scipy.signal import fftconvolve

        acf = fftconvolve(c_centered, c_centered[::-1], mode="full")
    except Exception:
        acf = np.correlate(c_centered, c_centered, mode="full")
    acf = acf[len(acf) // 2 :]  # 正延迟部分，acf[0] = sum of squares
    acf = acf / acf[0]  # 归一化 → rho_0 = 1

    tau = -1.0
    prev_pair = None
    half = len(acf) // 2
    for t in range(half):
        if 2 * t + 1 >= len(acf):
            break
        pair = acf[2 * t] + acf[2 * t + 1]
        if pair <= 0:
            break
        if prev_pair is not None and pair > prev_pair:
            break
        tau += 2.0 * pair
        prev_pair = pair

    tau = max(tau, 1e-12)
    return float(n / tau)
