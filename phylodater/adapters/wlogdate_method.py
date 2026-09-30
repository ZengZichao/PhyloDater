"""
WLogDateMethod - wLogDate 定年适配器

实现 wLogDate 定年功能。
wLogDate 通过最小化对数尺度突变率方差来定年系统发育树。

兼容性要点：
- CLI 参数：-i（输入树）、-t（校准文件）、-o（输出文件，必需）
- 向后时间模式（-b）：由 wLogDate 自行处理方向，PhyloDater 传入正向 Ma、
  读取正向 Ma（输出里的负值按冲突披露，不再取绝对值）
- 配置访问：适配器只读 WLogDateConfig 子配置（存于 ``self._wld_config``；
  ``self.config`` 则保留传入形态，与既有测试合约一致）
- 仅提供点估计（-p 参数控制优化初始点数量，非置信区间）
"""

import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ..core import DatingMethod, DatingMethodRegistry
from ..core.exceptions import (
    CalibrationError,
    ExecutionError,
    ResultParsingError,
    SemanticDegradationWarning,
)
from ..infrastructure import ProcessRunner, get_logger
from ..infrastructure.configuration import (
    CommonConfig,
    SoftwarePaths,
    ToolConfig,
    WLogDateConfig,
)
from ..infrastructure.safe_io import safe_writer
from ..models import (
    AgeConstraint,
    CalibrationPoint,
    CIType,
    DatingResult,
    NodeAgeEstimate,
    PhylogeneticTree,
)

# --------------------------------------------------------------------------- #
# wLogDate 输出树注解解析
# --------------------------------------------------------------------------- #
#
# wLogDate 输出形如 ``)Label[t=6.0, mu=0.03]:len``，其中 ', mu=' 的逗号后**有空格**。
# 早期正则只接受 ',mu='（无空格）变体，而单元测试使用的恰好是无空格样例，因此真实
# 上游输出（', mu='）一个注解都匹配不到，导致 parse_results 返回 0 个节点年龄（
# 集成测试在 benchmark 数据上复现）。这里统一放宽：
#   - 逗号后可跟任意空白（兼容 ',mu=' 与 ', mu='）；
#   - t= 之后到 ] 之间允许任意内容（兼容 mu= / CI= 以及字段顺序变化）。
_ANNOTATION_RE = re.compile(
    r"\)([^\[\]()]*?)\[t=(-?[\d.]+(?:[eE][-+]?\d+)?)"
    r"(?:[,\s][^\[\]]*?)?\](?=[;:,)]|$)"
)
# 从节点的原始注解块（如 ``[t=6.0, mu=0.03]``）里取回 t 值。
_AGE_IN_ANNOTATION_RE = re.compile(r"\[t=(-?[\d.]+(?:[eE][-+]?\d+)?)")


class _WLogDateTreeNode:
    """wLogDate 输出树的极简节点（自包含解析用，不依赖 ete3）。

    wLogDate 输出的 ``)Name[t=..., mu=...]`` 注解会让 ete3 的 Newick 解析器
    直接抛 ``'NoneType' object has no attribute 'add_child'``（在 Python 3.11+/
    ete3 3.1 上实测），因此 MRCA 回退改用本模块内的轻量解析器：只关心拓扑
    （父子、叶子名）与每个内部节点的原始注解块。
    """

    __slots__ = ("name", "annotation", "children", "parent")

    def __init__(self, name: str = "", annotation: Optional[str] = None) -> None:
        self.name = name
        self.annotation = annotation  # 原始注解块，如 "[t=6.0, mu=0.03]"；无则 None
        self.children: List["_WLogDateTreeNode"] = []
        self.parent: Optional["_WLogDateTreeNode"] = None


def _parse_wlogdate_newick(newick: str) -> _WLogDateTreeNode:
    """解析 wLogDate 输出树（带 ``)Name[t=age, mu=rate]:len`` 注解）。

    支持：未命名叶节点/内部节点、注解块、分支长度、根节点结尾 ';'。
    其余格式异常抛 :class:`ResultParsingError`，由调用方降级为仅标签匹配。
    """
    text = newick.strip()
    if text.endswith(";"):
        text = text[:-1]
    pos = 0
    n = len(text)

    def peek() -> str:
        return text[pos] if pos < n else ""

    def read_until(stops: str) -> str:
        nonlocal pos
        start = pos
        while pos < n and text[pos] not in stops:
            pos += 1
        return text[start:pos]

    def parse_clade(
        parent: Optional[_WLogDateTreeNode] = None,
    ) -> _WLogDateTreeNode:
        nonlocal pos
        node = _WLogDateTreeNode()
        node.parent = parent
        if peek() == "(":
            pos += 1
            while True:
                node.children.append(parse_clade(node))
                c = peek()
                if c == ",":
                    pos += 1
                    continue
                if c == ")":
                    pos += 1
                    break
                raise ResultParsingError(
                    f"wLogDate: unexpected token {c!r} in output tree"
                )
        node.name = read_until("[:,").strip()
        if peek() == "[":
            pos += 1
            ann = read_until("]")
            if peek() != "]":
                raise ResultParsingError("wLogDate: unbalanced '[' in output tree")
            pos += 1
            node.annotation = "[" + ann + "]"
        if peek() == ":":
            pos += 1
            read_until(",;)")  # 分支长度丢弃（本项目不需要）
        return node

    root = parse_clade()
    rest = text[pos:].strip()
    if rest:
        raise ResultParsingError(
            f"wLogDate: trailing characters {rest!r} after output tree"
        )
    return root


def _wlogdate_find_leaf(
    node: _WLogDateTreeNode, name: str
) -> Optional[_WLogDateTreeNode]:
    """按叶节点名查找节点（深度优先）。"""
    if not node.children and node.name == name:
        return node
    for child in node.children:
        found = _wlogdate_find_leaf(child, name)
        if found is not None:
            return found
    return None


def _wlogdate_lca(
    node_a: _WLogDateTreeNode, node_b: _WLogDateTreeNode
) -> Optional[_WLogDateTreeNode]:
    """两个节点的最近公共祖先（LCA）。"""
    ancestors = set()
    cur: Optional[_WLogDateTreeNode] = node_a
    while cur is not None:
        ancestors.add(cur)
        cur = cur.parent
    cur = node_b
    while cur is not None:
        if cur in ancestors:
            return cur
        cur = cur.parent
    return None


class WLogDateMethod(DatingMethod[Union[ToolConfig, WLogDateConfig]]):
    """
    wLogDate 定年适配器

    特性：
    - 仅提供点估计（无置信区间）
    - 使用所有有效叶节点（无 max_taxa 限制）
    - 统一时间单位处理
    - 支持向后时间模式（-b）
    - 派生参数（根锚点、探测到的比对长度）只留在适配器实例上，
      绝不写回共享的 config.wlogdate 对象
    """

    def __init__(
        self,
        config: Union[ToolConfig, WLogDateConfig],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        # wLogDate 是 7 个适配器里唯一把 ``self.config`` 保留为"传入形态"的：
        # 它从不读 ``self.config``（参数一律经 ``self._wld_config``），而
        # ``tests/unit/adapters/test_wlogdate_method.py::test_config_access_pattern``
        # 明断言 ``method.config is tool_config``。因此这里不归一
        # ``self.config``，只把"真正被读的那一份"单独命名为
        # ``_wld_config: WLogDateConfig``，静态类型仍然精确。
        # ``isinstance`` 与原先的 ``hasattr(config, "wlogdate")`` 对两种入参等价：
        # ``ToolConfig`` 必有 ``wlogdate``，``WLogDateConfig`` 必无。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            wld_config = config.wlogdate
            if common is None:
                common = config.common
        else:
            wld_config = config
        super().__init__(config, output_dir, software_paths, common_config=common)
        # 子配置引用（不是副本）：它"不写回共享配置对象"的约束靠的就是单一引用。
        self._wld_config: WLogDateConfig = wld_config

        self.logger = get_logger()
        self._input_files: Dict[str, Path] = {}
        self._calibrations: List[CalibrationPoint] = []
        self._warnings: List[Any] = (
            []
        )  # 收集需要写入 DatingResult/runtime_metadata 的警告
        # 实际生效的根锚点 / 比对长度：由 prepare_inputs 派生，仅存于本实例。
        # 不写回 self._wld_config（那是跨方法、跨调用共享的配置对象）。
        self._effective_root_time: Optional[float] = None
        self._effective_sequence_length: Optional[int] = None
        self._execution_seconds: float = 0.0

    @property
    def method_name(self) -> str:
        return "wlogdate"

    def validate_environment(self) -> bool:
        """检测 wLogDate 是否可用"""
        runner = ProcessRunner()
        custom_path = self.get_software_path("wlogdate_bin")
        executable = custom_path if custom_path else self._wld_config.wlogdate_bin

        available = runner.check_executable(executable)

        if available:
            version = runner.get_version(executable, "-v")
            # launch_wLogDate.py 没有 --version（-v 是 --verbose），因此以发行包
            # 元数据作为版本号来源；无法取得时只记录可执行文件本身。
            self.software_version = (
                version or self._package_version("wLogDate") or executable
            )
            self.logger.info(f"wLogDate detected: {version or executable}")
        else:
            self.logger.warning(f"wLogDate not found: {executable}")

        # launch_wLogDate.py 依赖 Python 包 logdate，仅存在可执行文件不够
        try:
            import logdate  # noqa: F401
        except ImportError:
            self.logger.warning(
                "wLogDate Python dependency 'logdate' is not installed. "
                "Install the 'wlogdate' distribution "
                "(https://pypi.org/project/wlogdate/, or bioconda 'wlogdate') "
                "or the method will fail at execution time."
            )

        return available

    def prepare_inputs(
        self,
        tree: PhylogeneticTree,
        calibrations: List[CalibrationPoint],
        alignment_path: Optional[Path] = None,
    ) -> Dict:
        """
        生成 wLogDate 输入文件

        1. 决定实际生效的根锚点与比对长度（只存本实例，不写回共享配置）
        2. 准备树文件（保留内部标签，单行无空行）
        3. 生成校准文件（LCA 格式：name=taxa1+taxa2\\t<age>）
        """
        self._calibrations = calibrations
        self._warnings = []

        # wLogDate 的随机初始化器在缺少根节点锚定时容易断言失败
        # （node.time > root_node.time）。如果校准配置中显式包含根节点，
        # 就把其年龄作为 -r 的候选锚点，提供稳定的初始化。
        #
        # 该派生值只写进适配器实例（self._effective_root_time），**不写回**
        # self._wld_config.root_time：那个对象是跨方法、跨调用共享的配置实例，
        # 就地改写会让用户显式设置的根锚点在第二次调用时仍然存在（B-28）。
        user_root_time = self._wld_config.root_time
        derived_root_time = self._derive_root_time(calibrations)
        self._effective_root_time = self._resolve_root_time(
            user_root_time, derived_root_time
        )

        # 获取序列长度（探测值同样只留在实例上，不写回配置）
        seq_len = self._wld_config.sequence_length or 1000
        if alignment_path and alignment_path.exists():
            from ..infrastructure import AlignmentMetadataExtractor

            extractor = AlignmentMetadataExtractor()
            metadata = extractor.extract(alignment_path)
            seq_len = metadata.sequence_length
            self.logger.info(
                f"wLogDate: detected alignment length {seq_len} (config value "
                f"{self._wld_config.sequence_length!r} is left unchanged)"
            )
        self._effective_sequence_length = seq_len

        # 准备树文件（保留内部标签）
        tree_file = self.work_dir / "input_tree.nwk"
        tree.write(tree_file)

        # wLogDate 的 launch_wLogDate.py 会逐行读取输入树文件并尝试解析每一行；
        # 如果文件末尾带有空行，会导致解析空字符串失败（UnexpectedEndOfStreamError）。
        # 因此这里把树文件整理为单行、无空行的内容。
        if tree_file.exists():
            raw_tree = tree_file.read_text()
            tree_lines = [
                line.strip() for line in raw_tree.splitlines() if line.strip()
            ]
            if tree_lines:
                tree_file.write_text(tree_lines[0])

        # 生成校准文件（LCA 格式：name=taxa1+taxa2\tage）
        constraints_file = self.work_dir / "constraints.txt"

        written = 0
        with safe_writer(constraints_file) as f:
            for cal in calibrations:
                # 优先 mrca_leaf_pair，回退 resolved_taxa（与 treePL/mcmctree 一致）
                taxa = (
                    list(cal.mrca_leaf_pair)
                    if cal.mrca_leaf_pair
                    else list(cal.resolved_taxa or [])
                )
                if len(taxa) < 2:
                    self.logger.warning(
                        f"Calibration '{cal.name}' has no resolved taxa/mrca_leaf_pair, skipping"
                    )
                    continue

                taxa_str = "+".join(taxa)
                constraint = cal.age_constraint

                if constraint is None:
                    self.logger.warning(
                        f"Calibration '{cal.name}' has no age constraint, skipping"
                    )
                    continue

                # 统一使用 to_software_format 生成 wLogDate 格式，再做单位转换
                # （直接传入校准点名与叶节点串，模型内部完成占位符渲染）。
                # 降级披露由 to_software_format 内部发出（wLogDate 只接受点估计），
                # 这里不再自行取中点，避免与集中实现漂移。
                formatted = constraint.to_software_format(
                    "wlogdate", name=cal.name, taxa_str=taxa_str
                )
                age = self._convert_time_unit(
                    self._age_from_wlogdate_line(formatted, cal.name)
                )
                f.write(f"{cal.name}={taxa_str}\t{age}\n")
                written += 1

        if written == 0:
            raise CalibrationError(
                "wLogDate: no valid calibrations written. "
                "Check that each calibration has mrca_leaf_pair or resolved_taxa "
                "with at least two taxa."
            )

        self._input_files = {"tree": tree_file, "constraints": constraints_file}

        self.logger.info(f"Generated wLogDate input: {tree_file}, {constraints_file}")
        return self._input_files

    def _derive_root_time(
        self, calibrations: List[CalibrationPoint]
    ) -> Optional[float]:
        """从根校准点推导 -r 锚点候选值（按配置的 time_unit 归一化）。

        校准年龄以 Ma 存储（fixed_age/max_age 均为 Ma），需按当前配置的 time_unit
        归一化，确保与用户显式配置的 root_time 采用同一单位约定（用户配置值即按
        time_unit 解释）。这样 execute 阶段的 _prepare_age（即 _convert_time_unit，
        仅做一次 Ga→Ma 转换）对任意 time_unit 都得到正确的 Ma 值，避免
        time_unit=='Ga' 时把已是 Ma 的校准年龄再 ×1000 造成 1000× 偏差。
        """
        for cal in calibrations:
            if cal.is_root_node and cal.age_constraint is not None:
                root_age = self._extract_age(cal.age_constraint, cal.name)
                if root_age is None:
                    continue
                tu = getattr(self._wld_config, "time_unit", "Ma")
                derived = root_age / 1000.0 if tu == "Ga" else root_age
                self.logger.info(
                    f"wLogDate: derived root anchor from calibration '{cal.name}': "
                    f"{derived} ({tu})"
                )
                return derived
        return None

    def _resolve_root_time(
        self, user_root_time: Optional[float], derived_root_time: Optional[float]
    ) -> Optional[float]:
        """决定实际生效的 -r 值：用户的显式设置优先，派生值只在缺失时补位。

        两者都在且不一致时，保留用户值并发出对照式警告（给出两个值），而不是
        静默用校准值顶掉用户锚点。
        """
        if derived_root_time is None:
            return user_root_time
        if user_root_time is None:
            return derived_root_time
        if abs(float(user_root_time) - float(derived_root_time)) <= 1e-9:
            self.logger.info(
                f"wLogDate: user root_time={user_root_time} agrees with the root "
                "calibration; using the user value."
            )
            return user_root_time

        msg = (
            f"wLogDate: user-set root_time={user_root_time} conflicts with the root "
            f"calibration age {derived_root_time}. Keeping the user value as the "
            "-r anchor; the root calibration is NOT applied. Fix the configuration "
            "if the calibration was meant to drive the root date."
        )
        self.logger.warning(msg, category=SemanticDegradationWarning)
        self._warnings.append(SemanticDegradationWarning(msg))
        return user_root_time

    @staticmethod
    def _age_from_wlogdate_line(formatted: str, cal_name: str) -> float:
        """从 to_software_format('wlogdate') 的输出里取回代表年龄。

        形如 ``Name=taxa1+taxa2\\t<age>``：年龄是最后一个制表符字段，因此按最右
        侧制表符切分（``rpartition``）而不是整体 split —— ``taxa_str`` 万一含制表
        符时 ``split("\\t")`` 会因段数 >2 抛 ValueError，把一个校准点名错误伪装成
        "格式不符"。渲染形状不符时显式报错而不是跳过。
        """
        head, _, age_str = formatted.rpartition("\t")
        if not head or not age_str.strip():
            raise CalibrationError(
                f"wLogDate: unexpected rendering of the age constraint for "
                f"calibration '{cal_name}': {formatted!r} (expected "
                f"'<name>=<taxa>\\t<age>')"
            )
        try:
            return float(age_str)
        except ValueError as e:
            raise CalibrationError(
                f"wLogDate: age constraint of calibration '{cal_name}' rendered to "
                f"a non-numeric value: {formatted!r}"
            ) from e

    def _extract_age(
        self,
        constraint: Optional[AgeConstraint],
        cal_name: str,
        taxa_str: str = "",
    ) -> Optional[float]:
        """从约束对象提取 wLogDate 使用的代表年龄（配置单位，通常 Ma）。

        生产调用点：根校准点 -> ``-r`` 锚点（见 _derive_root_time）。
        内部统一调用 to_software_format('wlogdate', name=..., taxa_str=...)，
        因此与校准文件走同一条降级路径（软界/区间/概率先验 -> 点估计 + 警告），
        不在适配器里重复实现降级规则。
        """
        if constraint is None:
            # 没有年龄界就没有可提取的代表年龄；调用方（_derive_root_time）本来就
            # 按 Optional 处理，这里返回 None 而不是在 None 上调渲染方法。
            return None
        formatted = constraint.to_software_format(
            "wlogdate", name=cal_name, taxa_str=taxa_str
        )
        return self._age_from_wlogdate_line(formatted, cal_name)

    def _convert_time_unit(self, age: float) -> float:
        """时间单位转换（Ga -> Ma）。

        ``time_unit`` 不是 ``WLogDateConfig`` 的声明字段（只有
        ``MCMCTreeConfig`` 有），因此必须按可选属性读取：缺失时按 "Ma"
        处理，等价于原来的 ``hasattr`` 守卫。
        """
        if str(getattr(self._wld_config, "time_unit", "Ma")) == "Ga":
            return age * 1000
        return age

    def _prepare_age(self, age: float) -> float:
        """单位转换（保留用于兼容旧调用；不再做方向取反）。

        wLogDate CLI 在 -b（向后时间）模式下会自行对输入年龄取反，因此
        PhyloDater 传入的校准年龄/根节点年龄都应保持为正向 Ma。
        """
        return self._convert_time_unit(age)

    def execute(self) -> bool:
        """执行 wLogDate 分析"""
        tree_file = self._input_files["tree"]
        constraints_file = self._input_files["constraints"]
        output_file = self.work_dir / "wlogdate_output.tre"

        # 确定可执行文件路径
        custom_path = self.get_software_path("wlogdate_bin")
        executable = custom_path if custom_path else self._wld_config.wlogdate_bin

        # 实际生效值：优先用 prepare_inputs 派生并留在实例上的值，未调用
        # prepare_inputs 时退回配置值（只读，不改配置）。
        seq_len = self._effective_sequence_length
        if seq_len is None:
            seq_len = self._wld_config.sequence_length
        root_time = self._effective_root_time
        if root_time is None:
            root_time = self._wld_config.root_time

        # 构建 wLogDate 命令
        cmd = [
            executable,
            "-i",
            str(tree_file),
            "-t",
            str(constraints_file),
            "-o",
            str(output_file),
            "-p",
            str(self._wld_config.num_replicates),
        ]

        # 算法参数：无条件显式传给上游。此前"仅在与 PhyloDater 内置默认值不同才传"
        # 的写法等于假定上游默认与本地默认一致；上游换版本/改默认时会静默拿到上游
        # 的新默认值，而配置里记录的仍是 PhyloDater 侧的假设值。
        if seq_len:
            cmd.extend(["-l", str(seq_len)])
        if self.common_config.seed is not None:
            cmd.extend(["-s", str(self.common_config.seed)])
        cmd.extend(["-m", str(self._wld_config.max_iter)])
        cmd.extend(["-u", str(self._wld_config.pseudocount)])
        cmd.extend(["-z", str(self._wld_config.zero_branch_len)])

        # 向后时间模式
        if self._wld_config.backward_time:
            cmd.append("-b")

        # 根节点和叶节点时间（保持正向 Ma，-b 模式由 wLogDate 自行取反）
        if root_time is not None:
            rt = self._prepare_age(root_time)
            cmd.extend(["-r", str(rt)])
        if self._wld_config.leaf_time is not None:
            lt = self._prepare_age(self._wld_config.leaf_time)
            cmd.extend(["-f", str(lt)])

        self.logger.info(f"wLogDate command: {' '.join(cmd)}")

        runner = ProcessRunner(
            cwd=self.work_dir, timeout=self.common_config.timeout or 7200
        )
        started = time.monotonic()
        try:
            result = runner.run(cmd)
        except FileNotFoundError as e:
            raise ExecutionError(
                f"wLogDate executable not found: {executable}. "
                f"Install the 'wlogdate' distribution "
                f"(https://pypi.org/project/wlogdate/, or bioconda "
                f"'wlogdate') before using this method."
            ) from e
        self._execution_seconds = time.monotonic() - started

        if result.returncode != 0:
            # wLogDate (fixed_init_lib.date_from_root_and_leaves) 在随机初始化时会
            # 选内部节点当"伪根"。如果选的节点是某个校准节点的后代，
            # 断言 node.time > root_node.time 会失败（wLogDate 库自身 bug）。
            # 这通常意味着校准与树拓扑/分支长度不兼容。
            err = (result.stderr or "") + (result.stdout or "")
            hint = ""
            if "ModuleNotFoundError" in err and "logdate" in err:
                hint = (
                    "The wLogDate launcher failed to import the 'logdate' Python package. "
                    "Install the 'wlogdate' distribution "
                    "(https://pypi.org/project/wlogdate/, or bioconda 'wlogdate') "
                    "before using the wlogdate method."
                )
            elif "assert" in err and "node.time" in err:
                hint = (
                    "This is a known limitation in wLogDate's initial date sampler "
                    "(not a PhyloDater bug): a randomly selected internal node was a "
                    "descendant of a calibration point, violating the node.time > "
                    "root.time assertion. The failure is random-seed dependent. You "
                    "can try a different --seed, add more calibration points, use "
                    "different calibrations, or choose another dating method."
                )
            raise ExecutionError(
                f"wLogDate execution failed (returncode={result.returncode}). "
                f"{hint}"
                f"\nstderr: {result.stderr[:800]}"
            )

        # 验证输出文件存在
        if not output_file.exists():
            raise ExecutionError(
                f"wLogDate output file not found: {output_file}. "
                f"stdout: {result.stdout[:500]}, stderr: {result.stderr[:500]}"
            )

        self._input_files["output_tree"] = output_file

        self.logger.success("wLogDate execution completed")
        return True

    def parse_results(self) -> DatingResult:
        """
        解析 wLogDate 输出

        wLogDate 输出格式：Newick 树，内部节点标签含 [t=X.XXX,mu=Y.YYY]。
        在 -b（向后时间）模式下，wLogDate 自身会把时间转为正向，
        因此 PhyloDater 解析时直接使用 t 值，无需再次取反。
        """
        output_tree = self._input_files.get("output_tree")

        if not output_tree or not output_tree.exists():
            raise ResultParsingError("wLogDate output tree not found")

        with open(output_tree, "r") as f:
            dated_tree_newick = f.read().strip()

        if not dated_tree_newick:
            raise ResultParsingError("wLogDate output tree is empty")

        # 解析节点年龄（包含未命名内部节点）
        node_ages = self._parse_tree_ages(dated_tree_newick)

        result = DatingResult(
            method_name=self.method_name,
            run_id=f"{self.method_name}_{self.work_dir.name}",
            dated_tree_newick=dated_tree_newick,
            node_ages=node_ages,
            raw_output_path=self.work_dir,
            execution_seconds=self._execution_seconds,
            metadata={
                "num_replicates": self._wld_config.num_replicates,
                "backward_time": self._wld_config.backward_time,
                # 实际生效值（而非配置里的假设值）：便于复现与事后核对
                "effective_sequence_length": self._effective_sequence_length,
                "effective_root_time": self._effective_root_time,
                "max_iter": self._wld_config.max_iter,
                "pseudocount": self._wld_config.pseudocount,
                "zero_branch_len": self._wld_config.zero_branch_len,
                "execution_time": self._execution_seconds,
            },
        )

        for warning in self._warnings:
            result.add_warning(warning)

        self.logger.success(f"Parsed wLogDate results: {len(node_ages)} node ages")
        return result

    def _resolve_parsed_age(self, age: float, node_desc: str) -> float:
        """校验从输出树读到的节点年龄；**不再对负值取绝对值**。

        本适配器的约定是：-b（向后时间）模式下 wLogDate 输出的 t 已是正向 Ma。
        在该约定下负值本不该出现，一旦出现即意味着方向、单位或校准↔拓扑冲突。
        把它折成正数等于把一个诊断信号伪装成一个看似合理的年龄（也是"祖先比
        后代年轻"最容易蒙混过关的形式），因此这里保留原值并记 ERROR + 结果级警告。
        """
        if age < 0:
            msg = (
                f"wLogDate reported a NEGATIVE age ({age}) for {node_desc}. Under the "
                "documented convention (output t is already forward-time Ma in -b "
                "mode) this signals a direction/unit/calibration-conflict problem. "
                "The value is kept AS IS (not abs()) so that downstream tables show "
                "it: check that the run really used -b, that input ages share one "
                "unit, and that the calibrations are compatible with the topology."
            )
            self.logger.error(msg)
            self._warnings.append(SemanticDegradationWarning(msg))
        return age

    def _parse_tree_ages(self, newick: str) -> Dict[str, NodeAgeEstimate]:
        """
        从 wLogDate 输出树解析年龄

        wLogDate 格式：)Label[t=age,mu=rate]:branch_len
        根节点：)Root[t=age,mu=rate];
        未命名内部节点：)[t=age,mu=rate]:branch_len

        wLogDate 在 -b（向后时间）模式下输出的 t 已经是正向 Ma，
        因此这里不做方向取反；负值按异常披露（见 _resolve_parsed_age）。
        """
        node_ages: Dict[str, NodeAgeEstimate] = {}

        # 匹配内部节点注释：)Label[t=age,mu=...]: 或 )Label[t=age,mu=...];
        # Label 可为空（未命名内部节点）。注意 wLogDate 上游输出 ', mu='（带空格），
        # 因此必须使用兼容空白的 _ANNOTATION_RE，而不是 ',mu=' 的窄版本。
        matches = list(_ANNOTATION_RE.finditer(newick))

        # 记录已使用的注释位置，避免被未命名节点重复提取
        used_positions: set = set()

        # 第一遍：按节点标签名精确匹配校准点（最可靠）
        cal_by_name = {cal.name: cal for cal in self._calibrations if cal.name}
        for match in matches:
            label = match.group(1).strip()
            if label and label in cal_by_name:
                age = self._resolve_parsed_age(
                    float(match.group(2)), f"calibration '{label}'"
                )
                node_ages[label] = NodeAgeEstimate(mean_age=age, ci_type=CIType.NONE)
                used_positions.add(match.start())
                del cal_by_name[label]

        # 第二遍：对仍未匹配的校准点，回退到 MRCA 拓扑匹配。
        # 不用 ete3：wLogDate 输出的 )Name[t=..., mu=...] 注解会让 ete3 的 Newick
        # 解析器直接抛 'NoneType' object has no attribute 'add_child'，因此这里使用
        # 模块内的轻量解析器（_parse_wlogdate_newick + _wlogdate_lca）。
        if cal_by_name:
            try:
                tree_root = _parse_wlogdate_newick(newick)
            except Exception as e:
                missing = ", ".join(sorted(cal_by_name))
                self.logger.warning(
                    f"wLogDate: could not parse the output tree for topological "
                    f"(MRCA) matching ({e}); calibrations without a matching node "
                    f"label will be missing from the results: {missing}."
                )
            else:
                for cal_name, cal in list(cal_by_name.items()):
                    if not cal.mrca_leaf_pair:
                        continue
                    tip1, tip2 = cal.mrca_leaf_pair
                    node1 = _wlogdate_find_leaf(tree_root, tip1)
                    node2 = _wlogdate_find_leaf(tree_root, tip2)
                    if node1 is None or node2 is None:
                        self.logger.warning(
                            f"wLogDate: leaf '{tip1}' or '{tip2}' not found in the "
                            f"output tree; cannot locate '{cal_name}' by MRCA."
                        )
                        continue
                    mrca = _wlogdate_lca(node1, node2)
                    if mrca is None or mrca.annotation is None:
                        self.logger.warning(
                            f"wLogDate: MRCA of '{tip1}'/'{tip2}' for '{cal_name}' "
                            f"carries no [t=...] annotation; skipping."
                        )
                        continue
                    ann_match = _AGE_IN_ANNOTATION_RE.search(mrca.annotation)
                    if not ann_match:
                        self.logger.warning(
                            f"wLogDate: no age found in the annotation "
                            f"{mrca.annotation!r} of the MRCA for '{cal_name}'."
                        )
                        continue
                    age = self._resolve_parsed_age(
                        float(ann_match.group(1)),
                        f"calibration '{cal_name}' (MRCA fallback)",
                    )
                    node_ages[cal_name] = NodeAgeEstimate(
                        mean_age=age, ci_type=CIType.NONE
                    )

        # 补充未命名内部节点
        self._add_unnamed_internal_nodes(newick, node_ages, used_positions)

        return node_ages

    def _add_unnamed_internal_nodes(
        self,
        newick: str,
        node_ages: Dict[str, NodeAgeEstimate],
        used_positions: Optional[set] = None,
    ) -> None:
        """补充未命名的内部节点

        wLogDate 输出树只为校准点标注名字和 [t=...,mu=...]；
        扫描所有尚未分配的节点注释，将其作为未命名内部节点加入结果。
        """
        used_positions = used_positions or set()

        idx = 1
        for match in _ANNOTATION_RE.finditer(newick):
            if match.start() in used_positions:
                continue

            label = match.group(1).strip()
            # 如果某个带标签节点尚未被匹配（MRCA 回退也可能没覆盖到），
            # 也视为未命名内部节点，避免遗漏。
            if label and label in node_ages:
                continue

            while f"internal_node_{idx}" in node_ages:
                idx += 1
            node_desc = f"unnamed internal node #{idx}"
            age = self._resolve_parsed_age(float(match.group(2)), node_desc)
            node_ages[f"internal_node_{idx}"] = NodeAgeEstimate(
                mean_age=age, ci_type=CIType.NONE
            )
            idx += 1


# 注册适配器
DatingMethodRegistry.register("wlogdate", WLogDateMethod)
