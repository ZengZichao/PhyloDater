"""
Tree Rooting - 根节点定根功能

支持多种定根方法：
- 中点定根 (Midpoint rooting)
- 外群定根 (Outgroup rooting)
- 最小方差定根
- MAD 定根（外部程序）

本模块**不导入 ete3**（审阅项 A-5/B-10）：ete3 3.1.x 导入标准库 ``cgi``，而
``cgi`` 已在 Python 3.13 被移除，因此 ete3 在 3.13+ 上根本无法导入。拓扑运算
（重挂、MRCA、根到叶距离）一律走 :class:`phylodater.models.tree.PhylogeneticTree`
自带的零依赖 Newick 游走器。

关于"树看起来已经有根"（审阅项 B-22）：Newick 里"根节点有 2 个子节点"只是**书写
约定**，不代表根放对了位置。用户显式请求某种定根方法时，就必须真的按该方法重根；
只有外群定根在"未提供外群"时才回退为原样返回（此时由调用方决定启发式策略）。
是否真的改变了根位置通过 :attr:`RootingResult.was_rerooted` 如实报告。
"""

import re
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import List, Optional, Tuple

from ..core.exceptions import TreeValidationError
from ..infrastructure.logging import get_logger
from ..models import PhylogeneticTree

#: 最内层（自身不含嵌套括号，但可以含逗号）的支树：``\([^()]*\)``
_NEWICK_INNERMOST_CLADE = re.compile(r"\([^()]*\)")
#: Newick/NHX 注释 ``[...]``（折叠形状前先剥掉）
_NEWICK_ANNOTATION = re.compile(r"\[[^\]]*\]")
#: 单个"标签[:枝长[:支持率]]"原子：标签**不得含空白**（日志句子会被这条挡掉），
#: 但允许 ``'quoted label with spaces'`` 这种带引号的写法；字符集只排除 Newick
#: 结构字符，好让中文名、带 ``%``/``/`` 的学名标签不被误判成"不是树"。
_NEWICK_NUMBER = r"\d+(?:\.\d*)?(?:[eE][-+]?\d+)?"
_NEWICK_ATOM = re.compile(r"(?:'[^']*'|[^\s,;():']+)?(?::" + _NEWICK_NUMBER + r"){0,2}")


def _looks_like_newick(text: str) -> bool:
    """判断一段文本"整体就是一棵 Newick 树"（审阅项 C-32 的形状闸门）。

    旧实现把 MAD 的整段 stdout 直接交给 ``PhylogeneticTree``，于是真正的成因
    （"MAD 没产出输出文件"）被推迟成一个毫不相干的解析错误。这里做 fullmatch 级
    的判据：只有一条语句（恰好一个分号、且在结尾）、括号配平、把最内层支树反复折叠
    成占位符后**整体**只剩"标签[:枝长]"的逗号串。折叠式化简可处理任意嵌套深度，
    比一次性正则可靠。
    """
    if not text:
        return False
    # 允许树被折行打印：先压缩空白，多行日志因此变成"带空格的原子"而被拒
    stripped = " ".join(text.split())
    if not stripped.endswith(";") or stripped.count(";") != 1:
        return False
    body = _NEWICK_ANNOTATION.sub("", stripped[:-1]).strip()
    if not body or body.count("(") != body.count(")"):
        return False
    while True:
        folded = _NEWICK_INNERMOST_CLADE.sub("N", body)
        if folded == body:
            break
        body = folded
    if "(" in body or ")" in body:  # 折叠不干净 = 括号结构本身不合法
        return False
    return all(_NEWICK_ATOM.fullmatch(atom) for atom in body.split(","))


class RootingMethod(Enum):
    """定根方法"""

    MIDPOINT = auto()  # 中点定根
    OUTGROUP = auto()  # 外群定根
    MINVAR = auto()  # 最小方差定根
    MAD = auto()  # 最小祖先偏差定根（依赖外部 MAD 程序）


@dataclass
class RootingResult:
    """定根结果"""

    rooted_tree: PhylogeneticTree
    method: RootingMethod
    root_edge_length: Optional[float] = None
    outgroup_taxa: Optional[List[str]] = None
    # B-22：让"实际上没有重根"这件事在结果对象里就看得见，而不是只留在日志里，
    # 这样可复现性报告能如实记录"用户请求了 midpoint、但根位置未改变"。
    was_rerooted: bool = True
    skipped_reason: Optional[str] = None


class TreeRootingService:
    """
    树定根服务

    提供多种定根方法
    """

    def __init__(self) -> None:
        self.logger = get_logger()

    def root_tree(
        self,
        tree: PhylogeneticTree,
        method: RootingMethod = RootingMethod.MIDPOINT,
        outgroup_taxa: Optional[List[str]] = None,
        output_path: Optional[Path] = None,
        mad_bin: str = "mad",
    ) -> RootingResult:
        """
        对树进行定根

        Args:
            tree: 输入树
            method: 定根方法
            outgroup_taxa: 外群分类单元（外群定根时使用）
            output_path: 输出文件路径
            mad_bin: MAD 可执行文件路径（MAD 定根时使用）
        """
        if method == RootingMethod.MIDPOINT:
            return self._midpoint_root(tree, output_path)
        elif method == RootingMethod.OUTGROUP:
            if not outgroup_taxa:
                raise ValueError("Outgroup taxa required for outgroup rooting")
            return self._outgroup_root(tree, outgroup_taxa, output_path)
        elif method == RootingMethod.MINVAR:
            return self._minvar_root(tree, output_path)
        elif method == RootingMethod.MAD:
            return self._mad_root(tree, mad_bin, output_path)
        else:
            raise ValueError(f"Unknown rooting method: {method}")

    # ------------------------------------------------------------------ #
    # 内部工具
    # ------------------------------------------------------------------ #

    def _as_model_tree(self, tree: PhylogeneticTree) -> PhylogeneticTree:
        """拿到一个确定带拓扑原语的 :class:`PhylogeneticTree`。

        调用方一般已经传的是真模型对象；对只带 ``.newick`` 的鸭子类型（旧的适配层
        调用点、测试替身）也照样能工作，而不是去依赖某个具体后端属性。
        """
        if isinstance(tree, PhylogeneticTree):
            return tree
        newick = getattr(tree, "newick", None)
        if isinstance(newick, str) and newick.strip():
            return PhylogeneticTree(newick)
        raise TreeValidationError(
            "定根服务需要一个带 Newick 文本的树对象（PhylogeneticTree 或任何有 "
            "非空 .newick 字符串的对象），实际收到 "
            f"{type(tree).__name__}。"
        )

    def _was_rerooted(
        self, original: PhylogeneticTree, rooted: PhylogeneticTree
    ) -> bool:
        """根位置是否真的改变：用有根同构指纹比较，而不是比 Newick 字符串。

        Newick 里两棵子树写在根的哪一侧没有含义，字符串相等不是可用判据；而
        "沿同一条枝把根往前挪 0.05"确实改变了树（各枝长随之重分配），指纹能抓住它。
        """
        try:
            return original._rooted_digest() != rooted._rooted_digest()
        except Exception:  # 判不出来时宁可如实说"变了"，也不谎报"没变"
            return True

    def _finish(
        self,
        original: PhylogeneticTree,
        rooted_newick_tree: Optional[PhylogeneticTree],
        method: RootingMethod,
        output_path: Optional[Path],
        message: str,
        outgroup_taxa: Optional[List[str]] = None,
        root_edge_length: Optional[float] = None,
        skipped_reason: Optional[str] = None,
    ) -> RootingResult:
        """统一收尾：判定是否真的改变了根位置、落盘、按结论记日志。

        ``rooted_newick_tree`` 为 ``None`` 表示计算失败/无需重根，此时原样返回输入对象
        （保持既有调用方对"返回同一个树对象"的依赖）。
        """
        rooted_tree = rooted_newick_tree if rooted_newick_tree is not None else original
        was_rerooted = rooted_newick_tree is not None and self._was_rerooted(
            original, rooted_tree
        )
        if rooted_newick_tree is None:
            reason = skipped_reason or "无法按请求的方法确定新的根位置"
            # 静默跳过正是 B-22 批评的形态：这里用 warning 而不是 info 级"看起来已有根"
            self.logger.warning(
                f"Tree was NOT re-rooted with {method.name}: {reason}. "
                f"输入树保持原样返回。"
            )
        elif not was_rerooted:
            self.logger.warning(
                f"{message}: 计算得到的根位置与输入树的现有根位置相同（"
                f"{skipped_reason or '该树已经根在这个分裂上'}）。"
            )
        else:
            self.logger.info(message)

        if rooted_newick_tree is not None and output_path:
            rooted_tree.write(output_path)

        return RootingResult(
            rooted_tree=rooted_tree,
            method=method,
            root_edge_length=root_edge_length,
            outgroup_taxa=outgroup_taxa,
            was_rerooted=was_rerooted,
            # 只要"根没动"，就必须带上理由——无论是算不出来，还是算出来恰好在此处。
            skipped_reason=(
                None
                if was_rerooted
                else (
                    skipped_reason
                    if rooted_newick_tree is None
                    else skipped_reason or "计算得到的根位置与现有根位置相同"
                )
            ),
        )

    # ------------------------------------------------------------------ #
    # 三种定根方法（零 ete3 依赖）
    # ------------------------------------------------------------------ #

    def _midpoint_root(
        self, tree: PhylogeneticTree, output_path: Optional[Path] = None
    ) -> RootingResult:
        """
        中点定根

        找到树中最长路径（tip-to-tip 距离最大的一对叶之间）的中点作为根节点。

        **方法学前提**：中点定根隐含"分子钟成立"——它把最长路径强行等分，等于在
        拟合速率之前先把"速率恒定"这个本应被检验的假设当成事实写进拓扑。因此
        :class:`phylodater.core.pipeline.DatingPipeline` 的默认 ``reroot_strategy``
        是 ``none``；有化石/外群信息时请优先用
        :attr:`RootingMethod.OUTGROUP`（或 MAD/MDCAT 等基于速率模型的定根）。

        B-22：只要用户显式请求本方法就**真的**重根。输入树以二叉形式书写
        （"看起来已有根"）不再构成跳过的理由——那只是 Newick 的书写约定。
        """
        model = self._as_model_tree(tree)
        try:
            rooted = model.rerooted_at_midpoint()
        except TreeValidationError:
            raise
        except Exception as e:  # 实现内部异常也必须响亮失败，不给"看似重了根"的树
            raise TreeValidationError(f"Midpoint rooting failed: {e}")

        if rooted is None:
            return self._finish(
                tree,
                None,
                RootingMethod.MIDPOINT,
                output_path,
                message="",
                skipped_reason=(
                    "该树的枝长全为 0/缺失，tip-to-tip 最长路径长度不超过 0，"
                    "中点位置在数学上不唯一"
                ),
            )

        return self._finish(
            tree,
            rooted,
            RootingMethod.MIDPOINT,
            output_path,
            message="Tree rooted at midpoint",
            skipped_reason="计算得到的中点就是现有的根位置",
        )

    def _outgroup_root(
        self,
        tree: PhylogeneticTree,
        outgroup_taxa: List[str],
        output_path: Optional[Path] = None,
    ) -> RootingResult:
        """
        外群定根

        以外群作为根节点的位置：新根放在"外群所在枝"上，外群整体作为根的两个子节点
        之一并保留其原枝长，因此所有 tip-to-tip 距离不变。多个外群取其 MRCA。

        关键修复（既有）：显式指定外群时必须尊重用户意图，即使当前树看似已定根
        （根有 2 个子节点）也要重新定根到指定外群；仅在未提供外群时才直接返回原树
        （此时由调用方决定启发式回退）。
        """
        if not outgroup_taxa:
            # 调用方（root_tree）已拦下；此处保留防御，语义与原实现一致。
            return self._finish(
                tree,
                None,
                RootingMethod.OUTGROUP,
                output_path,
                message="",
                outgroup_taxa=outgroup_taxa,
                skipped_reason="未提供外群分类单元",
            )

        model = self._as_model_tree(tree)
        known = set(model.tip_names)
        matched = [taxon for taxon in outgroup_taxa if str(taxon) in known]
        if not matched:
            raise TreeValidationError(
                f"None of the outgroup taxa found in tree: {outgroup_taxa}"
            )
        if len(matched) != len(outgroup_taxa):
            missing = [taxon for taxon in outgroup_taxa if str(taxon) not in known]
            self.logger.warning(
                f"部分外群分类单元不在树中，将只按已匹配到的 {matched} 定根；"
                f"未匹配：{missing}"
            )

        rooted, clade = model.rerooted_with_outgroup(matched)
        if rooted is None:
            return self._finish(
                tree,
                None,
                RootingMethod.OUTGROUP,
                output_path,
                message="",
                outgroup_taxa=outgroup_taxa,
                skipped_reason=(
                    f"外群 {sorted(clade) or matched} 的后代叶集合覆盖整棵树，"
                    "无法据此确定根位置（需要外群是真树的一个真子集）"
                ),
            )

        return self._finish(
            tree,
            rooted,
            RootingMethod.OUTGROUP,
            output_path,
            message=f"Tree rooted with outgroup: {outgroup_taxa}",
            outgroup_taxa=outgroup_taxa,
            skipped_reason="输入树本来就根在该外群所在的枝上",
        )

    def _minvar_root(
        self, tree: PhylogeneticTree, output_path: Optional[Path] = None
    ) -> RootingResult:
        """
        最小方差定根

        在树的现有节点中寻找"根到各叶距离方差最小"的位置定根（判据口径同
        ``numpy.var`` 的总体方差）。

        与中点定根同族：这一步同样隐含速率大致恒定，因此推荐优先使用外群定根。
        B-22：用户显式请求即真的重根，"看起来已有根"不构成跳过理由。
        """
        model = self._as_model_tree(tree)
        try:
            result = model.rerooted_at_min_variance()
        except TreeValidationError:
            raise
        except Exception as e:
            raise TreeValidationError(f"Minimum-variance rooting failed: {e}")

        if result is None:
            return self._finish(
                tree,
                None,
                RootingMethod.MINVAR,
                output_path,
                message="",
                skipped_reason="叶节点少于 2 个，根到叶距离方差无从定义",
            )

        rooted, variance = result
        return self._finish(
            tree,
            rooted,
            RootingMethod.MINVAR,
            output_path,
            message=(
                "Tree rooted at minimum variance position "
                f"(variance: {variance:.6f})"
            ),
            skipped_reason="方差最小的位置就是现有的根位置",
        )

    def _mad_root(
        self,
        tree: PhylogeneticTree,
        mad_bin: str = "mad",
        output_path: Optional[Path] = None,
    ) -> RootingResult:
        """
        使用 MAD (Minimal Ancestor Deviation) 定根

        Args:
            tree: 输入树
            mad_bin: MAD 可执行文件路径
            output_path: 输出文件路径

        Returns:
            定根结果

        输出定位（审阅项 C-32）：旧实现"猜一个 ``{stem}_rooted.nwk`` 文件名，猜不到
        就把整段 stdout 当 Newick"，于是"MAD 没写出文件"这一真实成因被推迟成遥远的
        解析错误。现在按三步**确定性**定位，三条失败路径各有独立错误信息：

        1. 把输出文件名显式告诉 MAD（第二个位置参数）；
        2. 运行被关在**本次专用的临时目录**里（``cwd`` 也是它），因此"目录里出现的
           别的文件"必然出自 MAD —— 即使该版本不认第二个参数、按自己的命名规则出文件，
           也能确定性地找到；
        3. 前两步都拿不到树，才回退 stdout，且回退前先过 :func:`_looks_like_newick`。
        """
        import shutil
        import subprocess
        import tempfile

        model = self._as_model_tree(tree)
        work_dir = Path(tempfile.mkdtemp(prefix="phylodater_mad_"))
        input_file = work_dir / "input.nwk"
        # 显式命名的输出文件：不再猜 ``{stem}_rooted.nwk``
        output_file = work_dir / "rooted.nwk"

        try:
            input_file.write_text(model.newick, encoding="utf-8")

            cmd = [mad_bin, str(input_file), str(output_file)]
            try:
                completed = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60,
                    cwd=str(work_dir),
                )
            except subprocess.TimeoutExpired:
                raise TreeValidationError("MAD timed out")

            if completed.returncode != 0:
                raise TreeValidationError(f"MAD failed: {completed.stderr}")

            new_newick, source = self._locate_mad_output(
                work_dir, input_file, output_file, completed.stdout
            )

            try:
                rooted_tree = PhylogeneticTree(new_newick)
            except Exception as e:
                # 与"没有输出文件"分开的另一条失败路径：拿到了形状正确的文本，但仍解析失败
                raise TreeValidationError(
                    f"MAD output from {source} looks like Newick but failed to parse: {e}"
                )

            # 定根是纯拓扑操作，不改叶集合；不一致说明读到的不是这棵树（例如误抓了
            # 别的产物文件），绝不能再往下当成"已定根的输入树"用。
            expected, observed = set(model.tip_names), set(rooted_tree.tip_names)
            if expected != observed:
                raise TreeValidationError(
                    f"MAD output from {source} is not the same tree as the input: "
                    f"tips {sorted(expected - observed)} disappeared, "
                    f"{sorted(observed - expected)} appeared."
                )

        except TreeValidationError:
            raise
        except Exception as e:
            raise TreeValidationError(f"Failed to reroot with MAD: {e}")
        finally:
            # 一次性临时目录整体删除：既不泄漏中间文件，也不会与并发运行的命名相撞
            shutil.rmtree(work_dir, ignore_errors=True)

        return self._finish(
            tree,
            rooted_tree,
            RootingMethod.MAD,
            output_path,
            message="Tree rooted with MAD",
            skipped_reason="MAD 给出的根位置与现有根位置相同",
        )

    def _locate_mad_output(
        self,
        work_dir: Path,
        input_file: Path,
        output_file: Path,
        stdout: str,
    ) -> Tuple[str, str]:
        """确定性定位 MAD 的定根树，返回 ``(Newick 文本, 来源说明)``。

        三条失败路径给三条不同的错误信息，其中"没有产出文件"直接点名真实成因。
        """
        if output_file.exists():
            text = output_file.read_text(encoding="utf-8", errors="replace").strip()
            if not _looks_like_newick(text):
                raise TreeValidationError(
                    f"MAD wrote the requested output file {output_file.name} but its "
                    f"content is not a Newick tree: {text[:160]!r}"
                )
            return text, output_file.name

        # 该 MAD 版本没按给定的名字出文件：在本次专用临时目录里找它留下的树
        # （目录私有 ⇒ 里面的文件必然来自 MAD）。候选按名字排序保证可复现，并把
        # 名字里带 "root" 的排在前面（多数版本的定根树就叫 *.rooted / *_rooted.nwk）；
        # 万一挑错文件，下面的"叶集合必须一致"会把问题响亮地报出来而不是静默接受。
        strays: List[Path] = []
        candidates = [
            p
            for p in work_dir.iterdir()
            if p != input_file and p.is_file() and p.name != output_file.name
        ]
        for candidate in sorted(
            candidates, key=lambda p: ("root" not in p.name.lower(), p.name)
        ):
            try:
                if candidate.stat().st_size == 0:
                    continue
                text = candidate.read_text(encoding="utf-8", errors="replace").strip()
            except OSError:
                continue
            if _looks_like_newick(text):
                return text, candidate.name
            strays.append(candidate)

        # 最后才回退 stdout——先过形状检查，别把日志喂给解析器
        text = (stdout or "").strip()
        if _looks_like_newick(text):
            return text, "stdout"

        raise TreeValidationError(
            "MAD exited successfully but produced no rooted tree: the requested "
            f"output file {output_file.name!r} was not created"
            + (
                f" and the {len(strays)} file(s) it did leave "
                f"({', '.join(p.name for p in strays)}) are not Newick"
                if strays
                else " and it wrote no other file"
            )
            + "; its stdout is not a Newick tree either "
            f"({(text[:160] + '...') if len(text) > 160 else (text or '<empty stdout>')!r}). "
            "This is a MAD invocation/output-convention problem, not a parsing bug — "
            "check the MAD binary and its expected output argument."
        )

    # ------------------------------------------------------------------ #
    # 只读查询（零 ete3 依赖）
    # ------------------------------------------------------------------ #

    def is_rooted(self, tree: PhylogeneticTree) -> bool:
        """检查树是否已按 Newick 约定定根（根节点恰有 2 个子节点）。

        注意这只是**书写形态**，不等于"根的位置是对的"（见 B-22）。
        解析失败时 :meth:`PhylogeneticTree.get_rooting_info` 保守返回未定根，
        因此本方法不会抛异常，也不会把"看不懂"说成"已定根"。

        C-31：旧实现在 ``ImportError``（ete3 装不上，本工作区的常态）时记一条 error
        然后 ``return True``，即"探测失败"被伪装成"已经定过根"这个最省事的答案，
        并与 B-22 的"看似已定根就跳过"叠加成静默跳过重根。现在探测链路零 ete3 依赖，
        且失败方向改为 fail-closed（答"未定根"，由上层按 ``require_rooted`` 处理）。
        """
        return bool(
            self._as_model_tree(tree).get_rooting_info().get("is_rooted", False)
        )

    def get_root_children(self, tree: PhylogeneticTree) -> List[List[str]]:
        """获取根节点的两个子分支各自的叶节点列表。

        返回 ``[左子树叶名列表, 右子树叶名列表]``；根节点不是二叉时返回 ``[]``
        （与原 ete3 实现的返回形态一致）。
        """
        model = self._as_model_tree(tree)
        nw = model._nw_tree
        if nw is None:
            return []
        children = nw.root.children
        if len(children) != 2:
            return []
        return [[label for label in child.leaf_labels() if label] for child in children]
