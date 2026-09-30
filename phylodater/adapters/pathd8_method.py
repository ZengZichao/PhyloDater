"""
PATHd8Method - PATHd8 定年适配器

实现 PATHd8 软件的适配器，支持固定年龄和区间约束
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Type, Union

from ..core import DatingMethod, DatingMethodRegistry
from ..core.exceptions import (
    CalibrationError,
    ExecutionError,
    PhyloDaterError,
    ResultParsingError,
    SemanticDegradationWarning,
)
from ..infrastructure import ProcessRunner, get_logger
from ..infrastructure.configuration import (
    CommonConfig,
    PATHd8Config,
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
    strip_leading_newick_comments,
)

# ---------------------------------------------------------------------------
# 上游 PATHd8 输出契约（下列串与列格式均已在随项目提供的
# `PhyloDater-参考软件/原始代码-PATHd8/`（含可执行二进制）上实测核对，
# 不是关键词猜测。见审阅报告 B-11 / B-12 / C-10。）
# ---------------------------------------------------------------------------

# PATHd8 的三类约束指令关键字（struct.h:201-203 R8S_FIX/R8S_MIN/R8S_MAX）
KIND_FIX = "fixage"
KIND_MIN = "minage"
KIND_MAX = "maxage"
_CONSTRAINT_KINDS = (KIND_FIX, KIND_MIN, KIND_MAX)

# NUM_FIX==0 时 meta.c:29-33 用 dprint(RES_FILE, ...) 打印到**输出文件**里：
# 它把用户的全部 min/max 约束清空（remove_all_constraints）并把根固定为 1.0
# 个任意时间单位。这是"结果单位已被改掉"的唯一硬证据。
NO_FIXNODE_MARKERS: Tuple[str, ...] = (
    "No fixnodes were defined",
    "Any user given constraint is ignored and root age is fixed to 1",
)

# 上游 error() → exit() 路径。ERR_FORMAT == 0（struct.h:36），也就是说
# `contradicting_node_error()` / `check_double_syntax()` 这类致命错误
# **以返回码 0 退出**，返回码不可判据（审阅报告 §六 末尾"共性风险"）。
# 这些串走 printf() → 标准输出，不在结果文件里，故必须同时扫描两者。
FATAL_OUTPUT_MARKERS: Tuple[str, ...] = (
    "This is a contradiction.",  # input.c:460, 482
    "Given constraints are contradicting each other.",  # d8.c:340
    "cannot be converted to a number",  # io_basics.c:124
    "Negative number",  # io_basics.c:128（"Negative number %f is not allowed"）
    "does not exist in the Newick tree",  # input.c:560
    'Exiting from function "',  # output.c:9
    "Wrong syntax",  # input.c:596
)

# 上游"我改写了你的约束"时唯一会说话的地方（input.c:448-475，printf → 标准输出）。
# 注意：真正对应"min/max 不相容 → update_constraints() 静默改写"的那句
# （meta.c:41）在上游被注释掉了，所以**关键词扫描无法覆盖该场景**，
# 必须改用下面的"约束对账"（informative 计数 + 节点表逐节点比对）。
UPSTREAM_REWRITE_MARKERS: Tuple[str, ...] = (
    "PATHd8 will consider the node to be",  # input.c:452/467/474
    "is defined both as a",  # input.c:472 fix_and_minmax_node_warning
    "is defined as a",  # input.c:466 double_minmax_warning
    "is defined to be",  # input.c:449/457/479 同节点多重定义
)

# output.c:246-248 —— 上游把"实际生效的约束数"写进结果文件（update_constraints /
# remove_all_constraints 都会把这些计数减下去），因此这是可用的对账依据。
_INFORMATIVE_COUNT_RE = re.compile(
    r"Number of informative\s+(fix|min|max)nodes\s*:\s*(-?\d+)", re.IGNORECASE
)

# 我们自己写出去的指令（§六.11 已逐字符核对语法正确，此处只按此语法回读）
_DIRECTIVE_RE = re.compile(
    r"^\s*mrca\s*:\s*(?P<t1>[^,;]+),\s*(?P<t2>[^,;]+),\s*"
    r"(?P<kind>fixage|minage|maxage)\s*=\s*(?P<value>[-+0-9.eE]+)\s*;",
    re.IGNORECASE | re.MULTILINE,
)
_NAME_DIRECTIVE_RE = re.compile(
    r"^\s*name of mrca\s*:\s*(?P<t1>[^,;]+),\s*(?P<t2>[^,;]+),\s*"
    r"name\s*=\s*(?P<name>[^;]+?)\s*;",
    re.IGNORECASE | re.MULTILINE,
)

# output.c:221 的节点表：
#   Ancestor of | Ancestor of | Name | Age(%10.3f) | #Terminals(%11d) |
#   MPL(%16.3f) | Rate(%14f) | minage(%11.1f 或 "-") | maxage(%11.1f 或 "-")
# "Age 三位小数 + 紧跟一个整数 + 再跟三位小数" 这一段足以把它和
# output.c:160 的 MPL 表（"MPL +/- 误差  终端数  Acc"）区分开。
_AGE_TABLE_ROW_RE = re.compile(
    r"^\s*\S+\s+\S+\s+(?P<name>\S+)\s+(?P<age>-?\d+\.\d{3})\s+(?P<nterm>\d+)\s+"
    r"(?P<mpl>-?\d+\.\d{3})\s+(?P<rate>\S+)\s+(?P<cmin>\S+)\s+(?P<cmax>\S+)\s*$",
    re.MULTILINE,
)
_AGE_TABLE_HEADER_RE = re.compile(
    r"^\s*Ancestor of\s+Ancestor of\s+Name\s+Age\s+#Terminals\s+MPL",
    re.IGNORECASE | re.MULTILINE,
)

# 根年龄落在该区间且用户从未给出 ≈1 Ma 的约束 → 疑似"根=1"的相对时间树
_UNNORMALIZED_ROOT_BAND = (0.5, 1.5)


class PATHd8Method(DatingMethod[PATHd8Config]):
    """
    PATHd8 定年适配器

    特性：
    - 要求至少一个固定年龄约束（并在**写盘之后**按产物复核，见 B-12）
    - 支持 minage/maxage 分开写入
    - 用"上游自报的生效约束计数 + 节点表"与**实际写出的输入文件**逐节点对账，
      检测 PATHd8 是否改写了用户约束（见 B-11；旧的 "conflict"/"adjusted"
      关键词检测在上游输出里逻辑上不可能命中，已删除）
    """

    def __init__(
        self,
        config: Union[ToolConfig, PATHd8Config],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        # 入口保留双形态（历史上两种调用都存在），但在交给基类之前就先归一：
        # 基类把 ``self.config`` 记为 ``PATHd8Config``，适配器后续只读子配置。
        # ``isinstance`` 与原先的 ``hasattr(config, "pathd8")`` 对这两种入参等价：
        # ``ToolConfig`` 必有 ``pathd8``，``PATHd8Config`` 必无。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            method_config = config.pathd8
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
        self._warnings: List[Any] = (
            []
        )  # 收集需要写入 DatingResult/runtime_metadata 的警告
        # B-12：约束对账以"实际写出去的文件"为准，这里缓存回读结果
        self._written_constraints: List[Dict[str, Any]] = []
        self._constraint_report: Dict[str, Any] = {}
        # B-11/C-10：上游的关键诊断信息走 printf → 标准输出，不在结果文件里，
        # 必须自己留住，否则 parse_results 无从判断这一轮是否真的跑通。
        self._process_stdout: str = ""
        self._process_stderr: str = ""
        self._sequence_length_meta: Dict[str, Any] = {}

    @property
    def method_name(self) -> str:
        return "pathd8"

    def _get_executable(self) -> str:
        """获取 PATHd8 可执行文件路径，优先使用 software_paths 自定义路径。"""
        path = self.get_software_path("pathd8_bin")
        return path if path else self.config.pathd8_bin

    def validate_environment(self) -> bool:
        """检测 PATHd8 是否可用"""
        from ..infrastructure.software_dependencies import SoftwareDependencyManager

        runner = ProcessRunner()
        executable = self._get_executable()
        available = runner.check_executable(executable)

        if available:
            version = runner.get_version(executable, "--version")
            # PATHd8 不支持 --version，也不在任何输出里写版本号（已核对上游
            # PATHd8.c / headers/output.c）。这里刻意**不**记录二进制的绝对安装
            # 路径：那既是机器专属信息（破坏跨机复现），又会把用户主目录一类的
            # 隐私路径泄漏到 scorecard.tsv / benchmark_results.json 等发布产物里。
            self.software_version = version or "PATHd8 (no version banner)"
            self.logger.info(f"PATHd8 detected: {version or 'available'}")
        else:
            self.logger.warning(f"PATHd8 not found: {executable}")
            install_cmd = SoftwareDependencyManager.get_install_command("pathd8")
            self.logger.warning(f"Installation suggestions:\n{install_cmd}")

        return available

    def _get_root_representative_leaves(
        self, tree: PhylogeneticTree
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
                    leaves.append(leaf.name)
            if len(leaves) == 2:
                return leaves[0], leaves[1]
        except Exception as e:
            self.logger.warning(f"Could not pick root representative leaves: {e}")

        names = tree.tip_names[:2]
        return names[0], names[1]

    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
        生成 PATHd8 输入文件

        1. 验证至少有一个固定年龄约束（"用户意图"层的快速失败）
        2. 安全化校准点名称
        3. 生成 infile
        4. **回读 infile**，按真正交给 PATHd8 的产物复核约束（B-12）

        第 4 步不可省：PATHd8 的 `fixage` 模式不是可选模式——`PATHd8.c:30-52`
        的 `main()` 无条件先跑 `run_mpl_fix()` 再跑 `run_mpl()`；而
        `meta.c:29-33` 在 `NUM_FIX == 0` 时会打印
        "Any user given constraint is ignored and root age is fixed to 1"、
        调用 `remove_all_constraints()` 清空用户的全部 min/max 约束，并把根
        固定为 1.0 个任意时间单位，**程序仍正常产出可解析的 `d8 tree`**
        （已用随项目提供的 PATHd8 二进制实测复现）。所以只看内存里的校准
        列表挡不住"唯一那个固定年龄校准因缺 mrca_leaf_pair 被跳过"这条路。
        """
        self._calibrations = calibrations
        self._written_constraints = []
        self._constraint_report = {}

        # 验证至少有一个固定年龄约束
        has_fixed = any(
            isinstance(cal.age_constraint, FixedAgeConstraint) for cal in calibrations
        )

        if not has_fixed:
            raise CalibrationError(
                "PATHd8 requires at least one fixed age (fixage) calibration. "
                "Please add a FixedAgeConstraint or use another software."
            )

        # 获取序列长度（C-12：不再凭空伪造默认值）
        seq_len, seq_len_source = self._resolve_sequence_length(alignment_path)
        self._sequence_length_meta = {
            "sequence_length": seq_len,
            "sequence_length_source": seq_len_source,
        }

        # 生成输入文件
        infile_path = self.work_dir / "PATHd8.infile"
        outfile_path = self.work_dir / "PATHd8.outfile"

        # 没有真正落盘的校准（B-6 的下游半段：上游加载器可能已经把点丢了）
        not_written: List[str] = []

        with safe_writer(infile_path) as f:
            # 序列长度（无比对时不写该指令，见 _resolve_sequence_length）
            if seq_len is not None:
                f.write(f"Sequence length = {seq_len};\n")

            # 树内容
            f.write(f"{tree.newick}\n\n")

            # PATHd8 不支持 is_root 关键字：用跨根的两个代表叶构造一个**覆盖根本身**
            # 的 mrca 约束。固定年龄写 fixage；区间/单边界约束同样可以落在根上
            # （已用 PATHd8 二进制实测：根 mrca 的 minage/maxage 会被计入
            # "informative minnodes/maxnodes"，根年龄被约束在该区间内）。旧实现只处
            # 理 FixedAgeConstraint，其余类型一律"不写文件、仅告警"，于是 ML 树上根
            # 年龄会被 MPL 推到约束区间之外（实测 153 Ma vs 请求的 100-110 Ma）。
            # 所有根校准统一在这里处理（C-11：过去 next(...) 只取第一个，
            # 而下面的循环又对所有 is_root_node 无条件 continue，
            # 于是第二个根校准既不写文件也不告警）。
            root_anchor_written = False
            for root_cal in [c for c in calibrations if c.is_root_node]:
                if root_anchor_written:
                    reason = (
                        "the root node has already been anchored by another "
                        "calibration"
                    )
                    msg = (
                        f"Root calibration '{root_cal.name}' was NOT written to the PATHd8 "
                        f"input file ({reason}); it has no effect on the result."
                    )
                    self.logger.warning(msg, category=SemanticDegradationWarning)
                    self._warnings.append(SemanticDegradationWarning(msg))
                    not_written.append(root_cal.name)
                    continue

                rt1, rt2 = self._get_root_representative_leaves(tree)
                safe_name = self._sanitize_name(root_cal.name)
                constraint = root_cal.age_constraint

                if constraint is None:
                    # 没有年龄界的根校准无法写成 fixage/minage/maxage：与下面
                    # “缺 mrca_leaf_pair 就跳过”同一条口径，记账后跳过，而不是
                    # 在 None 上调 to_software_format 碰 AttributeError。
                    msg = (
                        f"Root calibration '{root_cal.name}' was NOT written to the "
                        "PATHd8 input file (it has no age_constraint); it has no "
                        "effect on the result."
                    )
                    self.logger.warning(msg, category=SemanticDegradationWarning)
                    self._warnings.append(SemanticDegradationWarning(msg))
                    not_written.append(root_cal.name)
                    continue

                if isinstance(constraint, FixedAgeConstraint):
                    f.write(f"mrca: {rt1}, {rt2}, fixage={constraint.fixed_age};\n")
                else:
                    # 按统一的约束渲染入口写 min/max 指令（概率型约束在 constraints.py
                    # 内部已发出各自的降级警告）。
                    try:
                        formatted = constraint.to_software_format(
                            "pathd8", t1=rt1, t2=rt2
                        )
                    except Exception as exc:
                        formatted = ""
                        render_error = str(exc)
                    if not formatted:
                        reason = (
                            f"{type(constraint).__name__} could not be rendered as a "
                            f"PATHd8 mrca directive ({render_error})"
                        )
                        msg = (
                            f"Root calibration '{root_cal.name}' was NOT written to the "
                            f"PATHd8 input file ({reason}); it has no effect on the result."
                        )
                        self.logger.warning(msg, category=SemanticDegradationWarning)
                        self._warnings.append(SemanticDegradationWarning(msg))
                        not_written.append(root_cal.name)
                        continue
                    f.write(f"{formatted}\n")
                    if isinstance(constraint, SoftLowerBoundConstraint):
                        msg = (
                            f"Soft lower bound at the root degraded to a hard minage for "
                            f"PATHd8: {root_cal.name}"
                        )
                        self.logger.warning(msg, category=SemanticDegradationWarning)
                        self._warnings.append(SemanticDegradationWarning(msg))
                f.write(f"name of mrca: {rt1}, {rt2}, name={safe_name};\n")
                root_anchor_written = True

            # 校准约束
            for cal in calibrations:
                if cal.is_root_node:
                    continue
                if not cal.mrca_leaf_pair:
                    msg = f"Skipping calibration '{cal.name}' in PATHd8: no mrca_leaf_pair"
                    self.logger.warning(msg)
                    self._warnings.append(SemanticDegradationWarning(msg))
                    not_written.append(cal.name)
                    continue

                t1, t2 = cal.mrca_leaf_pair
                safe_name = self._sanitize_name(cal.name)

                constraint = cal.age_constraint
                if constraint is None:
                    msg = (
                        f"Skipping calibration '{cal.name}' in PATHd8: "
                        "no age_constraint"
                    )
                    self.logger.warning(msg)
                    self._warnings.append(SemanticDegradationWarning(msg))
                    not_written.append(cal.name)
                    continue

                # 统一使用约束对象的 to_software_format 生成 PATHd8 格式
                # （直接传入 MRCA 两叶名，模型内部完成占位符渲染）。
                formatted = constraint.to_software_format("pathd8", t1=t1, t2=t2)
                f.write(f"{formatted}\n")

                # 对于从概率分布降级而来的约束，to_software_format 内部已发出警告；
                # 但 PATHd8 需要保留原来的语义降级提示。这里统一为每种非原生支持的
                # 约束类型补充一条适配器级警告，便于用户追踪。
                if isinstance(constraint, SoftLowerBoundConstraint):
                    msg = f"Soft lower bound degraded to hard minage for PATHd8: {cal.name}"
                    self.logger.warning(msg, category=SemanticDegradationWarning)
                    self._warnings.append(SemanticDegradationWarning(msg))
                elif isinstance(constraint, MaximumAgeConstraint):
                    msg = f"Maximum age constraint degraded to maxage for PATHd8: {cal.name}"
                    self.logger.warning(msg, category=SemanticDegradationWarning)
                    self._warnings.append(SemanticDegradationWarning(msg))
                elif isinstance(
                    constraint,
                    (
                        SoftBoundsConstraint,
                        GammaPriorConstraint,
                        SkewNormalConstraint,
                        SkewTConstraint,
                    ),
                ):
                    msg = f"{type(constraint).__name__} degraded to uniform range for PATHd8: {cal.name}"
                    self.logger.warning(msg, category=SemanticDegradationWarning)
                    self._warnings.append(SemanticDegradationWarning(msg))

                # 命名 MRCA
                f.write(f"name of mrca: {t1}, {t2}, name={safe_name};\n")

        self._input_files = {"infile": infile_path, "outfile": outfile_path}

        # B-12：守卫从"意图层"下移到"产物层"
        self._guard_written_constraints(infile_path, calibrations, not_written)

        self.logger.info(f"Generated PATHd8 input file: {infile_path}")

        return self._input_files

    def _resolve_sequence_length(
        self, alignment_path: Optional[Path]
    ) -> Tuple[Optional[float], str]:
        """确定写入 infile 的 `Sequence length`（C-12：不再凭空造一个 1000）。

        上游 io_basics.c:142-153 的实现是::

            res = check_double_syntax(len);
            if (SEQUENCE_LENGTH_EXIST == TRUE) { res = (int)(res * SEQUENCE_LENGTH); }

        即该指令会把树上**每一条分支长度**乘以这个数并取整（缺省不乘），
        同时决定 `EDGE_LENS_ARE_INTEGERS`（是否做分子钟检验）。
        因此"用户没给比对就填 1000"等价于按一个臆造的比例重标输入树，
        短枝还可能被取整成 0；审阅报告 E-9 未定位到的用途现已定位。

        无比对时选择**不写该指令**：已用随项目提供的 PATHd8 二进制实测，
        缺省该指令时程序正常完成，且各节点年龄与写 `Sequence length = 500`
        时逐节点一致（MPL 定年对标度不变），只是 MPL 列按原始分支长度计。
        这件事会写进日志与结果元数据，不再悄悄伪造。
        """
        if alignment_path and Path(alignment_path).exists():
            from ..infrastructure import AlignmentMetadataExtractor

            extractor = AlignmentMetadataExtractor()
            metadata = extractor.extract(alignment_path)
            return metadata.sequence_length, "alignment"

        msg = (
            "No alignment available: the PATHd8 'Sequence length' directive is "
            "omitted, so branch lengths are used exactly as given in the tree "
            "(they are NOT rescaled by a fabricated default of 1000). "
            "Pass an alignment if your tree is in substitutions/site and you want "
            "PATHd8 to convert it to substitution counts."
        )
        self.logger.warning(msg, category=SemanticDegradationWarning)
        self._warnings.append(SemanticDegradationWarning(msg))
        return None, "not_provided_branch_lengths_used_as_given"

    def _guard_written_constraints(
        self,
        infile_path: Path,
        calibrations: List[CalibrationPoint],
        not_written: List[str],
    ) -> None:
        """回读真正写出的 infile，复核约束是否真的生效（B-12）。

        检查的是**文件内容**而不是内存里的 CalibrationPoint 列表，
        因此无论校准是在哪儿掉的（本适配器的跳过分支、models/tree.py 的
        B-6 静默丢弃、或 YAML 里就没写对），只要产物里没有可用的
        `fixage=` 指令就会在这里被拦住。
        """
        text = self._read_text(infile_path)
        entries = self._parse_written_constraints(text)
        self._written_constraints = entries

        counts = {kind: 0 for kind in _CONSTRAINT_KINDS}
        for entry in entries:
            counts[entry["kind"]] += 1

        self._constraint_report = {
            "calibrations_supplied": len(calibrations),
            "calibrations_not_written": list(not_written),
            "directives_written": dict(counts),
            "constrained_nodes_written": len(
                {(e["name"] or e["pair"]) for e in entries}
            ),
            "infile": str(infile_path),
        }

        if counts[KIND_FIX] == 0:
            msg = (
                f"PATHd8 input file '{infile_path.name}' contains NO 'fixage=' "
                f"directive ({len(calibrations)} calibration(s) requested, "
                f"{self._constraint_report['constrained_nodes_written']} effective). "
                "PATHd8's fixage stage runs on every invocation, and with zero fixnodes "
                'it prints "Any user given constraint is ignored and root age is '
                'fixed to 1", calls remove_all_constraints() and sets the root to '
                "1.0 arbitrary time units - i.e. every min/max calibration is void "
                "and the reported ages are relative, not Ma. Add at least one "
                "FixedAgeConstraint that resolves to an MRCA (mrca_leaf_pair), or "
                "use another software."
            )
            if not_written:
                msg += " Calibrations that never reached the input file: " + ", ".join(
                    sorted(set(not_written))
                )
            raise CalibrationError(msg)

        if counts[KIND_MIN] == 0 and counts[KIND_MAX] == 0:
            msg = (
                f"PATHd8 input file '{infile_path.name}' carries only "
                f"{counts[KIND_FIX]} fixage directive(s) and no minage/maxage at all. "
                "The dating will be anchored by fixed ages alone."
            )
            self.logger.warning(msg, category=SemanticDegradationWarning)
            self._warnings.append(SemanticDegradationWarning(msg))

        if len(calibrations) > self._constraint_report["constrained_nodes_written"]:
            msg = (
                f"PATHd8 reconciliation at input time: {len(calibrations)} "
                f"calibration(s) supplied, but only "
                f"{self._constraint_report['constrained_nodes_written']} node(s) "
                f"received a directive in the input file "
                f"(fixage={counts[KIND_FIX]}, minage={counts[KIND_MIN]}, "
                f"maxage={counts[KIND_MAX]})."
                + (
                    " Dropped: " + ", ".join(sorted(set(not_written)))
                    if not_written
                    else ""
                )
            )
            self.logger.warning(msg, category=SemanticDegradationWarning)
            self._warnings.append(SemanticDegradationWarning(msg))

    def _sanitize_name(self, name: str) -> str:
        """安全化名称，替换特殊字符"""
        # 替换所有非字母数字和下划线的字符
        safe = re.sub(r"[^A-Za-z0-9_]", "_", name)
        return safe[:50]  # 限制长度

    # ------------------------------------------------------------------
    # 约束对账：把"实际写出去的指令"和"PATHd8 自报的生效约束"逐节点比对
    # （B-11 的可用替代；B-12 的输出侧兜底；C-10 的致命信息）
    # ------------------------------------------------------------------
    @staticmethod
    def _read_text(path: Path) -> str:
        """读取文本文件（不存在/不可读时返回空串，由调用方决定如何失败）。"""
        try:
            with open(path, "r") as f:
                return f.read()
        except OSError:
            return ""

    @staticmethod
    def _mrca_pair_key(
        pair: Tuple[Optional[str], Optional[str]],
    ) -> Tuple[str, str]:
        """把一对叶名归一成排序后的 ``(t1, t2)`` 字典键。

        ``re.Match.group()`` 的标注是 ``str | Any``，直接 ``tuple(sorted(...))``
        只能得到变长元组 ``tuple[str | Any, ...]``，没法当
        ``Dict[Tuple[str, str], str]`` 的键用；这里显式拆成二元组。
        缺失的一侧按空串处理（上游正则不会匹配出 None，写在这里只是为了
        不把"不可能"当作类型上的保证）。
        """
        left = (pair[0] or "").strip()
        right = (pair[1] or "").strip()
        return (left, right) if left <= right else (right, left)

    def _parse_written_constraints(self, text: str) -> List[Dict[str, Any]]:
        """从 infile 文本里解析出真正写出去的约束指令。

        返回 [{"kind", "value", "pair", "name"}]，其中 `pair` 是排序后的
        (t1, t2)，`name` 是该 mrca 上最后一个 `name of mrca:` 标签
        （PATHd8 一个节点只保留一个名字，故取最后一个）。
        """
        names_by_pair: Dict[Tuple[str, str], str] = {}
        for m in _NAME_DIRECTIVE_RE.finditer(text):
            key = self._mrca_pair_key((m.group("t1"), m.group("t2")))
            names_by_pair[key] = m.group("name").strip()

        entries: List[Dict[str, Any]] = []
        for m in _DIRECTIVE_RE.finditer(text):
            kind = m.group("kind").lower()
            if kind not in _CONSTRAINT_KINDS:
                continue
            try:
                value = float(m.group("value"))
            except ValueError:
                continue
            pair = self._mrca_pair_key((m.group("t1"), m.group("t2")))
            entries.append(
                {
                    "kind": kind,
                    "value": value,
                    "pair": pair,
                    "name": names_by_pair.get(pair),
                }
            )
        return entries

    @staticmethod
    def _expected_bounds(
        entries: List[Dict[str, Any]],
    ) -> Dict[Any, Dict[str, Optional[float]]]:
        """按上游 `fill_fixage_data()` 的合并规则，算出每个节点应生效的上下界。

        同节点多条 minage 取最大、多条 maxage 取最小（input.c:385/424），
        fixage 同时钉住上下界（input.c:367-369），min==max 会被上游折成
        fixnode（input.c:391-397）——三种情形在这里都得到同一个 [lo, hi]，
        因此对账时不会把这些"上游会明说的合并"误判成改写。
        """
        grouped: Dict[Any, Dict[str, List[float]]] = {}
        for entry in entries:
            key = entry["name"] or entry["pair"]
            bucket = grouped.setdefault(key, {KIND_FIX: [], KIND_MIN: [], KIND_MAX: []})
            bucket[entry["kind"]].append(entry["value"])

        expected: Dict[Any, Dict[str, Optional[float]]] = {}
        for key, bucket in grouped.items():
            lo: Optional[float] = None
            hi: Optional[float] = None
            if bucket[KIND_FIX]:
                lo = hi = bucket[KIND_FIX][0]
            else:
                # 同节点多条同向边界：min 取最大、max 取最小（input.c:385/424）。
                # min == max 时上游会把它折成 fixnode 并 printf 报告（input.c:391-397），
                # 折完两列仍等于该值，所以这里无需特殊处理。
                if bucket[KIND_MIN]:
                    lo = max(bucket[KIND_MIN])
                if bucket[KIND_MAX]:
                    hi = min(bucket[KIND_MAX])
            expected[key] = {"min": lo, "max": hi}
        return expected

    @staticmethod
    def _reported_bounds_from_output(
        content: str,
    ) -> Tuple[Dict[str, Dict[str, Any]], Optional[Dict[str, Any]]]:
        """解析 output.c:221 的节点表（Name / Age / minage / maxage 列）。

        返回 `({name: {"age","min","max"}}, first_row)`。`"-"` 记为 None
        （= 上游声明该节点在这一列上没有约束）；同名节点以最后一次出现为准。
        `first_row` 是表里第一行数据，即根节点（`print_ancestor_age()` 先打
        根自己的行再递归子节点，output.c:216-227），用于判断结果是否还是
        "根=1"的相对时间标度。表解析不到时返回 ({}, None)。
        """
        rows: Dict[str, Dict[str, Any]] = {}
        first_row: Optional[Dict[str, Any]] = None
        header = _AGE_TABLE_HEADER_RE.search(content)
        if header is None:
            return rows, None
        for m in _AGE_TABLE_ROW_RE.finditer(content, header.end()):
            name = m.group("name")

            def _num(token: str) -> Optional[float]:
                if token in ("-", ""):
                    return None
                try:
                    return float(token)
                except ValueError:
                    return None

            row = {
                "name": name,
                "age": float(m.group("age")),
                "min": _num(m.group("cmin")),
                "max": _num(m.group("cmax")),
            }
            if first_row is None:
                first_row = row
            if name != "-":
                rows[name] = row
        return rows, first_row

    @staticmethod
    def _reported_informative_counts(content: str) -> Dict[str, int]:
        """解析 output.c:246-248 的 "Number of informative *nodes" 三行。"""
        reported = {KIND_FIX: 0, KIND_MIN: 0, KIND_MAX: 0}
        found = False
        for m in _INFORMATIVE_COUNT_RE.finditer(content):
            kind = {"fix": KIND_FIX, "min": KIND_MIN, "max": KIND_MAX}[
                m.group(1).lower()
            ]
            reported[kind] = int(m.group(2))
            found = True
        return reported if found else {}

    @staticmethod
    def _bounds_match(expected: Optional[float], reported: Optional[float]) -> bool:
        """比较一个边界值。节点表的 min/max 列用 %11.1f 打印，故容差取 0.06。

        方向是**单向**的：只检查"我要的边界有没有被丢掉或改掉"。
        上游额外给出一个我们没有要求的边界不算失效——那来自同一节点上的
        另一条校准（PATHd8 一个节点只保留一个名字）或 min==max 折成 fixnode
        的合并（input.c:391-397，且该情形上游自己会 printf 报告）。
        """
        if expected is None:
            return True
        if reported is None:
            return False
        return abs(expected - reported) <= max(0.06, 0.001 * abs(expected))

    def _reconcile_constraints(self, content: str) -> Dict[str, Any]:
        """把实际写出的约束与 PATHd8 自报的生效约束逐节点比对（B-11）。

        上游那句本该报告改写的 printf 被注释掉了（meta.c:36-43），
        所以"conflict"/"adjusted" 关键词扫描是逻辑上不可触发的死分支。
        这里改用两条**确实存在于真实输出里**的证据：
        1. `Number of informative fix/min/maxnodes`（被 update_constraints /
           remove_all_constraints 减计时会同步变小）；
        2. 节点表的 minage/maxage 两列，逐节点与本适配器写出去的边界对照。
        """
        entries = self._written_constraints
        report: Dict[str, Any] = dict(self._constraint_report)
        report["upstream_said_it_rewrote"] = []
        report["mismatches"] = []
        report["voided_nodes"] = []
        report["unlocated_nodes"] = []
        report["verified"] = False

        if not entries:
            report["note"] = (
                "no constraint directives were recovered from the PATHd8 input "
                "file; constraint reconciliation is unavailable for this run"
            )
            self._constraint_report = report
            return report

        combined = (
            content
            + "\n"
            + (self._process_stdout or "")
            + "\n"
            + (self._process_stderr or "")
        )
        for marker in UPSTREAM_REWRITE_MARKERS:
            if marker in combined:
                report["upstream_said_it_rewrote"].append(marker)

        reported_counts = self._reported_informative_counts(content)
        written_counts = {kind: 0 for kind in _CONSTRAINT_KINDS}
        for entry in entries:
            written_counts[entry["kind"]] += 1
        report["directives_written"] = dict(written_counts)
        report["informative_nodes_reported"] = dict(reported_counts)

        reported_bounds, root_row = self._reported_bounds_from_output(content)
        report["reported_root_age"] = root_row["age"] if root_row else None
        expected = self._expected_bounds(entries)
        report["nodes_written"] = len(expected)
        report["nodes_named_in_table"] = len(reported_bounds)

        voided: List[str] = []
        unlocated: List[str] = []
        if not reported_bounds:
            # 表没打印出来（PRINT_ANCESTOR 未开启等）→ 只能退回到计数对账，
            # 并明确说明逐节点核验不可用，而不是假装"没检测到问题"。
            report["note"] = (
                "PATHd8's per-node Age/minage/maxage table was not found in the "
                "output; only the informative-node counts could be reconciled"
            )
        else:
            report["verified"] = True
            for key, want in expected.items():
                name = key if isinstance(key, str) else None
                row = reported_bounds.get(name) if name else None
                if row is None:
                    unlocated.append(
                        f"{key}: written bound(s) "
                        f"[min={want['min']}, max={want['max']}] but no row of that "
                        "node name was found in PATHd8's table (PATHd8 keeps only one "
                        "name per node, so two calibrations sharing one MRCA collide)"
                    )
                    continue
                if not self._bounds_match(want["min"], row["min"]) or not (
                    self._bounds_match(want["max"], row["max"])
                ):
                    voided.append(
                        f"{key}: written bound(s) "
                        f"[min={want['min']}, max={want['max']}] but PATHd8 used "
                        f"[min={row['min']}, max={row['max']}]"
                    )
        report["voided_nodes"] = voided
        report["unlocated_nodes"] = unlocated
        report["mismatches"] = voided + unlocated

        if reported_counts:
            report["directive_count_deficit"] = sum(
                max(0, written_counts[k] - reported_counts.get(k, 0))
                for k in _CONSTRAINT_KINDS
            )

        self._constraint_report = report
        return report

    def _verify_constraints_were_honoured(self, content: str) -> List[Any]:
        """B-11 的可用替代 + B-12 的输出侧兜底。

        三道检查，全部基于真实存在于 PATHd8 输出里的东西：
        1. `No fixnodes were defined ... root age is fixed to 1`（meta.c:30，
           dprint 到结果文件）→ 结果单位已不是 Ma，直接中止。
        2. 根年龄落在 (0.5, 1.5) 且用户从没给过这个量级的校准
           → 疑似未归一化（B-12 修复建议 2）→ 中止。
        3. 逐节点对账：我们写出去的边界与节点表里的 minage/maxage 列不一致
           → 全部不一致 = 中止；部分不一致 = 响亮警告 + 差异清单入元数据。
        """
        findings: List[Any] = []

        def _fail(
            msg: str, error_cls: Type[PhyloDaterError] = CalibrationError
        ) -> None:
            self.logger.error(msg)
            raise error_cls(msg)

        # (1) 上游自报"全部约束作废、根固定为 1"
        hits = [marker for marker in NO_FIXNODE_MARKERS if marker in content]
        if hits:
            _fail(
                "PATHd8 reported "
                + " / ".join(repr(h) for h in hits)
                + ". With zero fixnodes PATHd8 calls remove_all_constraints() and "
                "fixes the root at 1.0 arbitrary time units, so every calibration "
                "you supplied has been ignored and the ages in this output are "
                "RELATIVE (root = 1), not Ma. The adapter refuses to return them. "
                "Give at least one calibration a FixedAgeConstraint that resolves "
                "to an MRCA (mrca_leaf_pair), or use another software."
            )

        report = self._reconcile_constraints(content)
        root_age = report.get("reported_root_age")
        written_values = [e["value"] for e in self._written_constraints]
        lo, hi = _UNNORMALIZED_ROOT_BAND
        near_one = [v for v in written_values if lo <= v <= hi]
        if root_age is not None and lo <= root_age <= hi and not near_one:
            _fail(
                f"PATHd8's dated tree has root age {root_age} while the smallest "
                f"calibration you supplied is "
                f"{min(written_values) if written_values else 'n/a'} Ma: the result "
                "looks like an un-normalized relative timescale (root fixed to 1), "
                "not absolute ages. Refusing to return it as Ma."
            )

        voided = report.get("voided_nodes") or []
        unlocated = report.get("unlocated_nodes") or []
        nodes_written = report.get("nodes_written") or 0

        if not report.get("verified"):
            msg = (
                "PATHd8 constraint reconciliation could NOT be performed for this "
                f"run ({report.get('note', 'no per-node table available')}). "
                "The adapter therefore cannot promise that the returned ages honour "
                "your calibrations: PATHd8 drops unsound min/max constraints "
                "silently (meta.c:36-43 has that printf commented out), so treat "
                "the result as unverified."
            )
            self.logger.warning(msg, category=SemanticDegradationWarning)
            findings.append(SemanticDegradationWarning(msg))
        elif voided and len(voided) >= nodes_written:
            _fail(
                "PATHd8 ignored or rewrote EVERY constraint this adapter wrote to "
                "its input file. Differences (written -> used):\n  "
                + "\n  ".join(voided)
                + "\nUpstream's own warning about unsound min/max data is commented "
                "out (meta.c:36-43), so the only way to catch this is comparing the "
                "per-node table with the input file. The returned ages do not "
                "reflect your calibrations."
            )
        elif voided or unlocated or report.get("upstream_said_it_rewrote"):
            details = "\n  ".join(voided + unlocated)
            msg = (
                "PATHd8 did not carry all constraints through unchanged "
                f"({len(voided)} node(s) with an altered or dropped bound, "
                f"{len(unlocated)} node(s) not locatable in the output table"
                + (
                    "; upstream said so on stdout: "
                    + ", ".join(report["upstream_said_it_rewrote"])
                    if report.get("upstream_said_it_rewrote")
                    else ""
                )
                + "). This means your min/max set was mutually unsound or two "
                "calibrations share one MRCA node. Diff (written -> used):\n  "
                + details
                + "\nUse this result with caution; the full diff is recorded in "
                "metadata['constraint_reconciliation']."
            )
            self.logger.warning(msg, category=SemanticDegradationWarning)
            findings.append(SemanticDegradationWarning(msg))
        else:
            self.logger.info(
                f"PATHd8 constraint reconciliation OK: all {nodes_written} "
                "written node bound(s) are reflected in the output table"
            )

        return findings

    def _restore_internal_node_names(self, newick: str) -> str:
        """将 PATHd8 输出树中的内部节点名恢复为原始校准名称，并确保根节点有标签。

        PATHd8 会把校准名中的非字母数字字符替换为下划线（如 SURF-12 → SURF_12），
        导致输出树节点名与输入配置不一致。本方法用原始校准名称重新标注内部 MRCA
        节点，使输出树与输入配置保持一致；同时补充根节点标签（PATHd8 默认不写）。
        """
        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)

            # 建立 sanitized_name → original_name 映射
            sanitized_to_original = {}
            root_name = None
            for cal in self._calibrations:
                if cal.name:
                    sanitized_to_original[self._sanitize_name(cal.name)] = cal.name
                    if cal.is_root_node:
                        root_name = cal.name

            # 遍历内部节点，若其 label 的 sanitized 形式命中某校准名，则恢复
            for node in tree.traverse("postorder"):
                if node.is_leaf():
                    continue
                original = sanitized_to_original.get(
                    self._sanitize_name(node.name or "")
                )
                if original:
                    node.name = original

            # ETE3 write(format=1) 会丢弃根节点名，需要手动注入
            newick_out: str = tree.write(format=1)
            if root_name:
                newick_out = self._ensure_root_label_in_newick(newick_out, root_name)
            return newick_out
        except ImportError:
            self.logger.warning(
                "ete3 not available, cannot restore PATHd8 internal node names"
            )
            return newick
        except Exception as e:
            self.logger.warning(f"Failed to restore PATHd8 internal node names: {e}")
            return newick

    def _ensure_root_label_in_newick(self, newick: str, root_name: str) -> str:
        """在 Newick 字符串末尾注入根节点标签（ETE3 format=1 会丢弃根节点名）。"""
        if not root_name or not newick:
            return newick
        newick = newick.strip()
        if not newick.endswith(";"):
            return newick
        if re.search(rf"\){re.escape(root_name)}(?:[\d.eE+-]*)?;\s*$", newick):
            return newick
        body = newick[:-1].rstrip()
        if body.endswith(")"):
            return f"{body}{root_name};"
        match = re.search(r"\)(:[\d.eE+-]+)$", body)
        if match:
            prefix = body[: match.start() + 1]
            suffix = match.group(1)
            return f"{prefix}{root_name}{suffix};"
        return newick

    def execute(self) -> bool:
        """执行 PATHd8 分析"""
        infile = self._input_files["infile"]
        outfile = self._input_files["outfile"]

        cmd = [self._get_executable(), str(infile.name), str(outfile.name)]

        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 3600
        )
        result = runner.run(cmd, expected_returncodes={0, 1})

        # 上游的致命诊断走 printf → 标准输出，不进结果文件；不留下来就永远看不见。
        self._process_stdout = result.stdout or ""
        self._process_stderr = result.stderr or ""

        # PATHd8 返回码说明（已用随项目提供的二进制实测核对）：
        # - 一次完全正常、约束全部生效的运行返回 **1**（"Calculation finished."）。
        # - 上游的 error() 用 ERR_* 宏当退出码，而 `ERR_FORMAT == 0`
        #   （struct.h:36），即"约束彼此矛盾/数值无法解析"这类**致命**错误
        #   以 **0** 退出（审阅报告 §六 末尾"共性风险"）。
        # 因此返回码本身既不报错也不报成功：真正的判据是
        # "结果文件里有没有可用的 d8 tree" + parse_results 里的约束对账。
        if result.returncode not in [0, 1]:
            raise ExecutionError(
                f"PATHd8 execution failed (returncode={result.returncode}): "
                f"{result.stderr}\n{result.stdout}"
            )

        # 返回码不在 {0,1} 之外，但 0/1 本身都不代表结果可信：先按上游真实
        # 打印的致命信息判一次（C-10：不再"发现疑似错误却不中止"）。
        self._check_fatal_upstream_output()

        if result.returncode == 1:
            msg = (
                "PATHd8 exited with code 1. This is the code the bundled build "
                "prints after an ordinary successful run ('Calculation finished.'), "
                "so it is neither an error nor a guarantee that your constraints "
                "survived; the constraint reconciliation in parse_results decides "
                "whether the ages honour the calibrations you supplied."
            )
            self.logger.warning(msg)
            self._warnings.append(SemanticDegradationWarning(msg))

        self.logger.success("PATHd8 execution completed")
        return True

    def _check_fatal_upstream_output(self) -> None:
        """上游一声明致命就中止（C-10：原来的子串猜测 fail-open 已删）。

        旧代码是 `if "cannot be converted" in content or "Error" in content:`
        后只 `logger.warning("Potential error detected…")` 就继续解析：
        大小写敏感的 "Error" 既漏（上游打的是小写 `error("fill_fixage_data",ERR_FORMAT)`）
        又误伤正文，而且"发现疑似错误却不中止"使坏结果一路到底。
        现在按上游实际语句（FATAL_OUTPUT_MARKERS）匹配，命中即抛错。
        """
        combined = (self._process_stdout or "") + "\n" + (self._process_stderr or "")
        hits = [marker for marker in FATAL_OUTPUT_MARKERS if marker in combined]
        if not hits:
            return
        line = ""
        for candidate in combined.splitlines():
            if any(marker in candidate for marker in hits):
                line = candidate.strip()
                break
        raise ExecutionError(
            "PATHd8 reported a fatal input/format problem and aborted "
            f"(markers: {', '.join(hits)}). Note that PATHd8 exits with code 0 for "
            "ERR_FORMAT, so the exit code cannot be used as a success criterion. "
            f"Upstream message: {line!r}"
        )

    def parse_results(self) -> DatingResult:
        """解析 PATHd8 输出"""
        # 上游在 stdout 上的致命声明优先于一切（返回码不可判据：ERR_FORMAT == 0，
        # 见 execute 的注释）。放在结果文件检查之前，这样用户看到的是
        # PATHd8 自己的诊断，而不是"文件不存在"。
        self._check_fatal_upstream_output()

        outfile = self._input_files["outfile"]

        if not outfile.exists():
            raise ResultParsingError(f"PATHd8 output file not found: {outfile}")

        content = self._read_text(outfile)
        if not content.strip():
            raise ResultParsingError(
                f"PATHd8 output file is empty: {outfile}. PATHd8 aborts several fatal "
                "errors with exit code 0, so an empty result file is how a failed "
                "run looks here."
            )

        # 从 infile 回读"实际写出去的约束"（B-12 的产物层口径；prepare_inputs
        # 已缓存一份，这里允许 parse_results 被单独调用时重新读取）
        if not self._written_constraints:
            infile = self._input_files.get("infile")
            if infile is not None:
                self._written_constraints = self._parse_written_constraints(
                    self._read_text(Path(infile))
                )

        # 提取 d8 tree
        match = re.search(r"^d8 tree\s*:\s*(.+)$", content, re.MULTILINE)
        if not match:
            raise ResultParsingError("Could not find 'd8 tree' in PATHd8 output")

        dated_tree_newick = match.group(1).strip()

        # B-11/B-12：与上游自报的生效约束逐节点对账（取代过去那段
        # "conflict"/"adjusted" 关键词扫描——上游对应的 printf 被注释掉了，
        # 那个分支在真实输出上逻辑上不可能触发）。
        reconciliation_warnings = self._verify_constraints_were_honoured(content)

        # 将 PATHd8 输出树中的内部节点名恢复为原始校准名称
        # PATHd8 会把名称中的 '-' 等字符替换为 '_'，导致输出树节点名与输入配置不一致
        dated_tree_newick = self._restore_internal_node_names(dated_tree_newick)

        # 提取分子钟检验结果
        # PATHd8 输出包含极具生物学价值的信息：
        # - 对每个节点进行的分子钟检验 (Clock test) 结果
        # - 接受 (Accepted) 和拒绝 (Rejected) 的统计
        # - 每个节点的平均路径长度 (MPL)
        metadata = {}

        # 提取总体统计
        clock_tests_match = re.search(
            r"Clock tests\s*:\s*(\d+)\s*Accepted\s*:\s*(\d+)\s*Rejected\s*:\s*(\d+)",
            content,
        )
        if clock_tests_match:
            metadata["clock_tests"] = {
                "total": int(clock_tests_match.group(1)),
                "accepted": int(clock_tests_match.group(2)),
                "rejected": int(clock_tests_match.group(3)),
                "acceptance_rate": (
                    int(clock_tests_match.group(2)) / int(clock_tests_match.group(1))
                    if int(clock_tests_match.group(1)) > 0
                    else 0.0
                ),
            }

            # 记录到日志，帮助用户评估数据集的速率异质性
            acceptance_rate = metadata["clock_tests"]["acceptance_rate"]
            if acceptance_rate < 0.5:
                self.logger.warning(
                    f"PATHd8 clock test acceptance rate is low ({acceptance_rate:.1%}). "
                    f"This dataset may deviate significantly from the strict clock model."
                )
            else:
                self.logger.info(
                    f"PATHd8 clock test: {metadata['clock_tests']['accepted']}/"
                    f"{metadata['clock_tests']['total']} accepted ({acceptance_rate:.1%})"
                )

        # 提取每个节点的 MPL (Mean Path Length) 信息
        # 格式示例: "NodeName MPL=0.1234"
        mpl_pattern = r"(\w+)\s+MPL[=:]\s*([\d.]+)"
        mpl_matches = re.findall(mpl_pattern, content)
        if mpl_matches:
            metadata["mpl_values"] = {name: float(mpl) for name, mpl in mpl_matches}

        # 可复现性：写出去的 Sequence length 是哪来的、约束对账结果如何（C-12/B-11）
        metadata.update(self._sequence_length_meta)
        metadata["constraint_reconciliation"] = self._constraint_report

        # 解析节点年龄
        node_ages = {}

        # 从输出中提取各校准点的年龄
        for cal in self._calibrations:
            # 查找该校准点的年龄
            # PATHd8 输出格式示例: "LUCA 3500.00"
            pattern = rf"{re.escape(self._sanitize_name(cal.name))}\s+([\d.]+)"
            age_match = re.search(pattern, content)

            if age_match:
                age = float(age_match.group(1))
                node_ages[cal.name] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)

        # 如果无法从输出解析，使用树结构计算
        if not node_ages:
            node_ages = self._extract_ages_from_tree(dated_tree_newick)
        else:
            # PATHd8 的输出通常不包含根节点年龄，但树结构可以计算出来。
            # 如果存在根节点校准且尚未提取，则补上一个根年龄。
            root_cal = next(
                (
                    cal
                    for cal in self._calibrations
                    if cal.is_root_node and cal.name not in node_ages
                ),
                None,
            )
            if root_cal is not None:
                try:
                    from ete3 import Tree

                    tree = Tree(
                        strip_leading_newick_comments(dated_tree_newick), format=1
                    )
                    root = tree.get_tree_root()
                    root_age = max(
                        root.get_distance(leaf) for leaf in tree.get_leaves()
                    )
                    node_ages[root_cal.name] = NodeAgeEstimate(
                        mean_age=root_age, ci_type=CIType.NONE
                    )
                except Exception as e:
                    self.logger.warning(
                        f"Could not extract root age for {root_cal.name}: {e}"
                    )

        # 补充未命名的内部节点
        self._add_unnamed_internal_nodes(dated_tree_newick, node_ages)

        result = DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=0.0,  # 将在外层填充
            metadata=metadata,
        )

        # 添加警告
        for warning in reconciliation_warnings:
            result.add_warning(warning)
        for warning in self._warnings:
            result.add_warning(warning)

        self.logger.success(f"Parsed PATHd8 results: {len(node_ages)} node ages")

        return result

    def _add_unnamed_internal_nodes(
        self, newick: str, node_ages: Dict[str, NodeAgeEstimate]
    ) -> None:
        """补充未命名的内部节点

        PATHd8 输出的 d8 tree 通常只为校准点标注内部节点名；
        遍历整棵树，把其余内部节点的年龄也提取出来，使报告能显示全部内部节点。
        """
        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)
            root = tree.get_tree_root()
            root_age = max(root.get_distance(leaf) for leaf in tree.get_leaves())

            calibrated_names = {self._sanitize_name(name) for name in node_ages.keys()}
            idx = 1
            for node in tree.traverse("postorder"):
                if node.is_leaf():
                    continue
                node_name = self._sanitize_name(node.name) if node.name else ""
                # 已命名的校准节点（含根节点）不再重复加入
                if node_name and node_name in calibrated_names:
                    continue
                if node is root and "Root" in node_ages:
                    continue
                age = root_age - root.get_distance(node)
                if age < 0:
                    continue
                while f"internal_node_{idx}" in node_ages:
                    idx += 1
                node_ages[f"internal_node_{idx}"] = NodeAgeEstimate(
                    mean_age=age, ci_type=CIType.NONE
                )
                idx += 1
        except Exception as e:
            self.logger.warning(
                f"Could not extract unnamed internal nodes for PATHd8: {e}"
            )

    def _extract_ages_from_tree(self, newick: str) -> Dict[str, NodeAgeEstimate]:
        """从带分支长度的树中提取节点年龄

        PATHd8 输出的 d8 tree 的分支长度代表绝对时间
        """
        node_ages = {}

        try:
            from ete3 import Tree

            tree = Tree(strip_leading_newick_comments(newick), format=1)
            root = tree.get_tree_root()
            root_age = max(root.get_distance(leaf) for leaf in tree.get_leaves())

            for cal in self._calibrations:
                # 根节点校准：直接使用根年龄
                if cal.is_root_node:
                    node_ages[cal.name] = NodeAgeEstimate(
                        mean_age=root_age, ci_type=CIType.NONE
                    )
                    continue

                if not cal.mrca_leaf_pair:
                    self.logger.warning(
                        f"Skipping calibration '{cal.name}' in PATHd8 result parsing: no mrca_leaf_pair"
                    )
                    continue

                tip1, tip2 = cal.mrca_leaf_pair

                try:
                    mrca = tree.get_common_ancestor(tip1, tip2)
                    # PATHd8 树中分支长度为绝对时间
                    # root.get_distance(mrca) = 从根到MRCA的路径长度 = root_age - age
                    # age = 根年龄 - 根到MRCA距离
                    root_to_mrca = root.get_distance(mrca)
                    age = root_age - root_to_mrca

                    if age < 0:
                        self.logger.warning(
                            f"Negative age computed for {cal.name}: "
                            f"root_age={root_age}, root_to_mrca={root_to_mrca}"
                        )
                        continue

                    node_ages[cal.name] = NodeAgeEstimate(
                        mean_age=age, ci_type=CIType.NONE
                    )
                except Exception as e:
                    self.logger.warning(f"Could not extract age for {cal.name}: {e}")

            # 补充未命名的内部节点
            self._add_unnamed_internal_nodes(newick, node_ages)

        except ImportError:
            self.logger.warning("ete3 not available, cannot extract ages from tree")

        return node_ages

    def cleanup(self, preserve_intermediates: bool = False) -> Dict[str, Path]:
        """清理临时文件"""
        return {}


# 注册适配器
DatingMethodRegistry.register("pathd8", PATHd8Method)
