"""
InputValidator - 输入文件验证服务

支持多种主流格式的树文件和序列文件的存在性验证和格式有效性检查。
"""

import re
from dataclasses import dataclass, field
from enum import Enum
from io import StringIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from ..infrastructure.logging import get_logger

#: 格式判据的两种形态：``bytes`` 前缀，或已编译的 ``bytes`` 正则（按首行匹配）。
Signature = Union[bytes, "re.Pattern[bytes]"]


class TreeFormat(Enum):
    """支持的树文件格式"""

    NEWICK = "newick"
    NEXUS = "nexus"
    PHYLOXML = "phyloxml"
    STOCKHOLDT = "stockholdt"  # Stockholm 格式（如 RAxML 输出）
    AUTO = "auto"


class AlignmentFormat(Enum):
    """支持的序列比对文件格式"""

    FASTA = "fasta"
    PHYLIP = "phylip"
    PHYLIP_INTERLEAVED = "phylip_interleaved"
    CLUSTAL = "clustal"
    STOCKHOLM = "stockholm"
    AUTO = "auto"


@dataclass
class ValidationResult:
    """验证结果"""

    is_valid: bool
    format_detected: Optional[str] = None
    # ``= None`` 作为 ``List[str]`` 的默认值是类型上的谎言（运行期靠调用方
    # 自己防止 None 参拼）；``default_factory`` 才是那个真实行为。
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    # 中性信息（如"包含 N 个叶节点/序列"），不应以 WARNING 级别打印，
    # 避免给用户造成"出了问题"的误导。
    info: List[str] = field(default_factory=list)
    tree_count: Optional[int] = None  # 树文件中的树数量

    def __post_init__(self) -> None:
        if self.errors is None:
            self.errors = []
        if self.warnings is None:
            self.warnings = []
        if self.info is None:
            self.info = []


class InputValidator:
    """
    输入文件验证服务

    支持多种主流格式的树文件和序列文件的存在性验证和格式有效性检查。
    """

    # 树文件格式的魔术字节/特征（由 :meth:`_detect_tree_format` 唯一读取）
    #
    # 判据表必须显式标注：不标时 mypy 会把 ``list[bytes]`` 与
    # ``list[Pattern[bytes]]`` 归并成 ``object``，于是下面的迭代和传参全部落空。
    TREE_FORMAT_SIGNATURES: Dict[TreeFormat, List[Signature]] = {
        TreeFormat.NEXUS: [
            b"#NEXUS",
            b"#nexus",
        ],
        TreeFormat.PHYLOXML: [
            b"<phyloxml",
            b"<?xml",
        ],
    }

    # 序列文件格式的魔术字节/特征
    #
    # 审阅项 C-46 附带项：这张表此前是**死数据**——``_detect_alignment_format``
    # 一处硬编码 ``b">"`` 判 FASTA、另一处又把 PHYLIP 头部正则内联重编了一遍，
    # 于是表里的 FASTA / PHYLIP 两项从不被读取，只会误导维护者"改这里就够了"。
    # 现在它是比对格式判据的**唯一来源**（bytes 按前缀匹配，已编译的正则按整行
    # 匹配），改判据只需要动这张表。
    ALIGNMENT_FORMAT_SIGNATURES: Dict[AlignmentFormat, List[Signature]] = {
        AlignmentFormat.FASTA: [
            b">",
        ],
        AlignmentFormat.PHYLIP: [
            re.compile(rb"^\s*\d+\s+\d+\s*$"),  # PHYLIP 头部: "  num_taxa  seq_length"
        ],
        AlignmentFormat.CLUSTAL: [
            b"CLUSTAL",
            b"clustal",
        ],
        AlignmentFormat.STOCKHOLM: [
            b"# STOCKHOLM",
            b"#stockholm",
        ],
    }

    # 树文件扩展名 → 格式（与上面的签名表同理，是 C-46"死数据会误导维护者"的
    # 同族问题：这两张表此前**从不被读取**，检测函数里各自又硬编码了一遍扩展名
    # 清单，两处会漂移。现在扩展名判定只读这两张表。）
    TREE_FORMAT_BY_EXTENSION = {
        ".nwk": TreeFormat.NEWICK,
        ".newick": TreeFormat.NEWICK,
        ".tree": TreeFormat.NEWICK,
        ".treefile": TreeFormat.NEWICK,
        ".tre": TreeFormat.NEWICK,
        ".nxs": TreeFormat.NEXUS,
        ".nexus": TreeFormat.NEXUS,
        ".phyloxml": TreeFormat.PHYLOXML,
        ".xml": TreeFormat.PHYLOXML,
        # RAxML 的 .sth 信息文件不是树；检测出来之后由 _validate_tree_format
        # 明确拒绝（见审阅项 C-46），而不是悄悄当作 Newick 放行
        ".sth": TreeFormat.STOCKHOLDT,
        ".stockholm": TreeFormat.STOCKHOLDT,
    }

    ALIGNMENT_FORMAT_BY_EXTENSION = {
        ".fa": AlignmentFormat.FASTA,
        ".fasta": AlignmentFormat.FASTA,
        ".fna": AlignmentFormat.FASTA,
        ".faa": AlignmentFormat.FASTA,
        ".fas": AlignmentFormat.FASTA,
        ".afa": AlignmentFormat.FASTA,  # aligned fasta
        ".phy": AlignmentFormat.PHYLIP,
        ".phylip": AlignmentFormat.PHYLIP,
        ".phylip-interleaved": AlignmentFormat.PHYLIP_INTERLEAVED,
        ".aln": AlignmentFormat.CLUSTAL,
        ".clustal": AlignmentFormat.CLUSTAL,
        ".sto": AlignmentFormat.STOCKHOLM,
        ".stockholm": AlignmentFormat.STOCKHOLM,
    }

    # 向后兼容：扩展名清单由上面的表**派生**，不再是一份会各自漂移的独立清单
    TREE_EXTENSIONS = frozenset(TREE_FORMAT_BY_EXTENSION)
    ALIGNMENT_EXTENSIONS = frozenset(ALIGNMENT_FORMAT_BY_EXTENSION)

    def __init__(self) -> None:
        self.logger = get_logger()

    def validate_tree_file(
        self, file_path: Path, expected_format: TreeFormat = TreeFormat.AUTO
    ) -> ValidationResult:
        """
        验证树文件

        Args:
            file_path: 树文件路径
            expected_format: 期望的格式（AUTO 表示自动检测）

        Returns:
            ValidationResult: 验证结果
        """
        # 基础检查
        base_result = self._validate_file_exists(file_path)
        if not base_result.is_valid:
            return base_result

        # 格式检测
        if expected_format == TreeFormat.AUTO:
            detected_format = self._detect_tree_format(file_path)
            if detected_format is None:
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"无法自动检测树文件格式: {file_path}",
                        "请确保文件是有效的 Newick、Nexus 或 PhyloXML 格式",
                    ],
                )
        else:
            detected_format = expected_format

        # 格式验证
        format_result = self._validate_tree_format(file_path, detected_format)
        if not format_result.is_valid:
            return format_result

        # 尝试解析
        parse_result = self._try_parse_tree(file_path, detected_format)
        if not parse_result.is_valid:
            return parse_result

        return ValidationResult(
            is_valid=True,
            format_detected=detected_format.value,
            warnings=base_result.warnings
            + format_result.warnings
            + parse_result.warnings,
            info=base_result.info + format_result.info + parse_result.info,
        )

    def validate_alignment_file(
        self, file_path: Path, expected_format: AlignmentFormat = AlignmentFormat.AUTO
    ) -> ValidationResult:
        """
        验证序列比对文件

        Args:
            file_path: 序列比对文件路径
            expected_format: 期望的格式（AUTO 表示自动检测）

        Returns:
            ValidationResult: 验证结果
        """
        # 基础检查
        base_result = self._validate_file_exists(file_path)
        if not base_result.is_valid:
            return base_result

        # 格式检测
        if expected_format == AlignmentFormat.AUTO:
            detected_format = self._detect_alignment_format(file_path)
            if detected_format is None:
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"无法自动检测序列比对文件格式: {file_path}",
                        "请确保文件是有效的 FASTA、PHYLIP 或 Clustal 格式",
                    ],
                )
        else:
            detected_format = expected_format

        # 格式验证
        format_result = self._validate_alignment_format(file_path, detected_format)
        if not format_result.is_valid:
            return format_result

        # 尝试解析
        parse_result = self._try_parse_alignment(file_path, detected_format)
        if not parse_result.is_valid:
            return parse_result

        return ValidationResult(
            is_valid=True,
            format_detected=detected_format.value,
            warnings=base_result.warnings
            + format_result.warnings
            + parse_result.warnings,
            # info 此前被丢掉（树那条路是合并的），C-47 要求"报告实际使用的编码"
            # 必须能抵达调用方，这里对齐树文件的处理
            info=base_result.info + format_result.info + parse_result.info,
        )

    def validate_calibration_file(self, file_path: Path) -> ValidationResult:
        """
        验证校准配置文件

        审阅项 B-24：此前这道门只做五步结构检查（存在 / 扩展名 / YAML 可解析 /
        顶层是 dict / 有 ``calibrations`` 列表）就返回 ``is_valid=True``，CLI 随即
        打印"校准配置文件格式有效"。实测五份最普通的误写 YAML（``type`` 放错层级、
        ``min``/``max`` 放错层级、条目写成标量、空条目、没有 ``constraint`` 键）
        **全部通过**。现在这里做的是**内容级**校验，与
        :mod:`phylodater.services.calibration_loader` 共用同一份"合法类型 × 必填键"
        契约（避免两处漂移）：

        * 条目必须是映射且非空，``name`` 必须是非空字符串；
        * ``constraint`` 必须存在、必须是映射、必须有合法 ``type``；
        * 每类约束的必填键齐全、未知键（通常是拼错/放错层级）报错；
        * 所有年龄值必须是**有限**数值且 ``0 <= age <= 4600``，``min < max``；
          ``.nan`` / ``.inf`` 这类 YAML 字面量在解析层就被拒（同时收掉 B-19 的入口）；
        * 报错文案带条目序号与文件路径。

        Args:
            file_path: 校准配置文件路径

        Returns:
            ValidationResult: 验证结果
        """
        # 基础检查
        base_result = self._validate_file_exists(file_path)
        if not base_result.is_valid:
            return base_result

        # 检查文件扩展名
        ext = file_path.suffix.lower()
        if ext not in [".yaml", ".yml"]:
            return ValidationResult(
                is_valid=False,
                errors=[f"校准配置文件必须是 YAML 格式 (.yaml 或 .yml): {file_path}"],
            )

        # 尝试解析 YAML
        try:
            import yaml

            with open(file_path, "r", encoding="utf-8", newline="") as f:
                data = yaml.safe_load(f)

        except yaml.YAMLError as e:
            return ValidationResult(is_valid=False, errors=[f"YAML 解析错误: {e}"])
        except (OSError, UnicodeDecodeError) as e:
            return ValidationResult(
                is_valid=False, errors=[f"读取校准配置文件失败: {e}"]
            )

        # 结构 + 内容级校验（共用加载器里的同一份契约表）
        from .calibration_loader import (
            hard_errors,
            informationals,
            validate_calibration_payload,
        )

        messages = validate_calibration_payload(data, source=str(file_path))
        errors = hard_errors(messages)
        notes = informationals(messages)

        if errors:
            self.logger.error(
                f"校准配置文件内容校验失败 ({file_path}): {len(errors)} 处问题"
            )
            return ValidationResult(
                is_valid=False,
                format_detected="yaml",
                errors=errors,
                warnings=base_result.warnings + notes,
            )

        entries = (data or {}).get("calibrations") or []
        info = [f"校准配置文件包含 {len(entries)} 个校准条目"]
        return ValidationResult(
            is_valid=True,
            format_detected="yaml",
            warnings=base_result.warnings + notes,
            info=info,
        )

    def validate_taxonomy_file(self, file_path: Path) -> ValidationResult:
        """
        验证分类学信息文件

        支持格式：
        1. 两列格式（名称 + 分类学字符串）
        2. 多列 TSV 格式（name, domain, phylum, ...）

        审阅项 C-39：这里除了"首行是否两列 / 第二列是否以 ``d__`` 开头"的浅检查，
        还会跑一次 :meth:`DeepValidator.validate_taxonomy_file`（**分类学名称在
        不同层级间循环引用**的检测，如 ``d__A;p__B`` 与 ``d__B;p__A`` 并存）。
        那个校验器带完整的 ``TaxonomyConflictError`` 与可执行建议，但此前在生产
        路径上**零调用点**——CLI 只调本方法，于是"承重件已造好、线没接"。

        Args:
            file_path: 分类学信息文件路径

        Returns:
            ValidationResult: 验证结果
        """
        # 基础检查
        base_result = self._validate_file_exists(file_path)
        if not base_result.is_valid:
            return base_result

        # 尝试读取文件
        try:
            with open(file_path, "r", encoding="utf-8", newline="") as f:
                lines = []
                for i, line in enumerate(f):
                    lines.append(line.strip())
                    if i >= 10:  # 只读取前11行用于验证
                        break

            if not lines:
                return ValidationResult(
                    is_valid=False, errors=[f"分类学信息文件为空: {file_path}"]
                )

            # 检测分隔符
            first_line = lines[0]
            delimiter = "\t" if "\t" in first_line else ","

            # 检查列数
            parts = first_line.split(delimiter)
            if len(parts) < 2:
                return ValidationResult(
                    is_valid=False,
                    errors=[f"分类学信息文件格式错误: 至少需要两列: {file_path}"],
                )

            # 检查是否是两列格式（分类学字符串）
            if len(parts) == 2:
                if not self._has_taxonomy_string(lines, delimiter):
                    return ValidationResult(
                        is_valid=False,
                        errors=[
                            f"分类学信息文件格式错误: 第二列应以 'd__' 或 'd_' 开头"
                            f"（允许首行为表头）: {file_path}"
                        ],
                    )
            else:
                # 多列格式：检查是否有必要的列
                header = lines[0].lower()
                required_cols = ["name", "domain"]
                for col in required_cols:
                    if col not in header:
                        return ValidationResult(
                            is_valid=False,
                            errors=[f"分类学信息文件缺少必要的列 '{col}': {file_path}"],
                        )

        except Exception as e:
            return ValidationResult(
                is_valid=False, errors=[f"读取分类学信息文件失败: {e}"]
            )

        # 审阅项 C-39：把深度校验（循环依赖检测）接进生产路径
        deep_errors, deep_warnings, deep_info = self._deep_validate_taxonomy(file_path)
        if deep_errors:
            self.logger.error(
                f"分类学信息文件深度校验失败 ({file_path}): {len(deep_errors)} 处问题"
            )
            return ValidationResult(
                is_valid=False,
                format_detected="taxonomy_table",
                errors=deep_errors,
                warnings=base_result.warnings + deep_warnings,
                info=deep_info,
            )

        return ValidationResult(
            is_valid=True,
            format_detected="taxonomy_table",
            warnings=base_result.warnings + deep_warnings,
            info=deep_info,
        )

    @staticmethod
    def _has_taxonomy_string(lines: List[str], delimiter: str = "") -> bool:
        """表格里是否真的存在 ``d__…`` 式分类学字符串列（允许首行是表头）。

        审阅项 C-39 的前置条件：旧的浅检查只看**首行第二列**，于是
        ``examples/taxonomy.tsv``（首行 ``name\\ttaxonomy`` 是表头）这类完全合法的
        两列表格会被判"格式错误"，接在后面的深度校验永远跑不到。现在 sampled
        的前若干行里只要有任意一行第二列以 ``d__`` / ``d_`` 开头即认为形状正确；
        仍然是保守判定（不接受纯垃圾文本），也不放宽到"猜列名"。

        分隔符逐行判定（与 :meth:`DeepValidator._parse_taxonomy_lines` 同一取向：
        先试制表符再试逗号），避免"表头用逗号、数据用制表符"时误判。
        """
        for line in lines:
            if not line or line.startswith("#"):
                continue
            row_delimiter = "\t" if "\t" in line else (delimiter or ",")
            parts = line.split(row_delimiter)
            if len(parts) < 2:
                continue
            second = parts[1].strip().strip('"')
            if second.startswith("d__") or second.startswith("d_"):
                return True
        return False

    def _deep_validate_taxonomy(
        self, file_path: Path
    ) -> Tuple[List[str], List[str], List[str]]:
        """跑一次 :meth:`DeepValidator.validate_taxonomy_file`（审阅项 C-39）。

        那个校验器会检测分类学名称在不同层级之间的**循环引用**
        （``d__A;p__B`` 与 ``d__B;p__A`` 并存）并以 :class:`TaxonomyConflictError`
        上抛，附带可执行建议——但报告实测它在生产路径上零调用点。

        本方法的职责是把这条线接上，同时守住验证门的契约：**返回错误清单而不是
        抛异常**。CLI 的 ``_validate_input_files`` 不 catch 异常，直接抛出去会让
        "校验"变成裸 traceback；翻译成语义等价的 ``errors`` 后，调用方（以及
        ``--taxonomy-file`` 的所有既有通路）只需按 ``is_valid`` 分支处理。

        Returns:
            (errors, warnings, info): 三个消息列表，errors 非空即判定为无效
        """
        from ..core.exceptions import TaxonomyConflictError
        from .deep_validator import DeepValidator

        errors: List[str] = []
        warnings: List[str] = []
        info: List[str] = []

        try:
            detail = DeepValidator().validate_taxonomy_file(file_path)
        except TaxonomyConflictError as e:
            # 循环依赖：CRITICAL，必须拦下
            errors.append(f"分类学信息文件深度校验失败: {e}")
            if getattr(e, "suggestion", None):
                errors.append(str(e.suggestion))
            return errors, warnings, info
        except Exception as e:  # pragma: no cover - 防御：深度校验自身崩了不该判死
            self.logger.warning(
                f"分类学信息文件深度校验未能完成 ({file_path}) "
                f"{type(e).__name__}: {e}；本轮仅结构检查生效"
            )
            warnings.append(
                f"分类学循环依赖检查未能完成（{type(e).__name__}: {e}），"
                f"本次仅做了表头/列数检查"
            )
            return errors, warnings, info

        errors.extend(detail.errors)
        warnings.extend(detail.warnings)
        if detail.circular_dependencies:
            # 校验器只在非空时抛异常；这里再兜一层，保证"拿到循环清单"必然判无效
            errors.append(
                "分类学信息文件存在循环依赖: "
                + ", ".join(f"'{a}' <-> '{b}'" for a, b in detail.circular_dependencies)
            )
        if not errors:
            info.append(
                f"循环依赖检查覆盖 {detail.entries_checked} 条分类学记录（{file_path}）"
            )
        return errors, warnings, info

    def _validate_file_exists(self, file_path: Path) -> ValidationResult:
        """
        验证文件是否存在、是否为空

        Args:
            file_path: 文件路径

        Returns:
            ValidationResult: 验证结果
        """
        warnings = []

        if not file_path.exists():
            return ValidationResult(is_valid=False, errors=[f"文件不存在: {file_path}"])

        if not file_path.is_file():
            return ValidationResult(
                is_valid=False, errors=[f"路径不是文件: {file_path}"]
            )

        file_size = file_path.stat().st_size
        if file_size == 0:
            return ValidationResult(
                is_valid=False, errors=[f"文件为空 (0 字节): {file_path}"]
            )

        if file_size < 10:
            warnings.append(f"文件非常小 ({file_size} 字节): {file_path}")

        return ValidationResult(is_valid=True, warnings=warnings)

    # 读文本时的编码尝试顺序（审阅项 C-47：旧代码的 warning 承诺"将尝试其他编码"
    # 却没有任何代码真的尝试，且最终把原始编解码报错丢给用户）
    TEXT_ENCODINGS = ("utf-8-sig", "utf-8", "cp1252", "latin-1")

    def _read_text_best_effort(
        self, file_path: Path, limit: Optional[int] = None
    ) -> Tuple[Optional[str], Optional[str]]:
        """按 :data:`TEXT_ENCODINGS` 依次尝试读取文本。

        Returns:
            (text, encoding) 或 (None, None)：所有编码都失败时。
            ``encoding`` 报告的是**实际生效**的编码：``utf-8-sig`` 编解码器在没有
            BOM 时与 ``utf-8`` 完全等价，因此只有文件头真的是
            ``b'\\xef\\xbb\\xbf'`` 才回报 ``utf-8-sig``，否则普通的 UTF-8 文件会被
            误报成"检测到 UTF-8 BOM"（C-47 要求披露真实编码，就不能虚构 BOM）。
        """
        for encoding in self.TEXT_ENCODINGS:
            try:
                with open(file_path, "r", encoding=encoding, newline="") as f:
                    text = f.read(limit) if limit else f.read()
                if encoding == "utf-8-sig" and not self._starts_with_utf8_bom(
                    file_path
                ):
                    encoding = "utf-8"
                return text, encoding
            except (UnicodeDecodeError, LookupError):
                continue
            except OSError:
                return None, None
        return None, None

    @staticmethod
    def _starts_with_utf8_bom(file_path: Path) -> bool:
        """文件头是否为 UTF-8 BOM（读失败按 False 处理，交给后续检查报错）。"""
        try:
            with open(file_path, "rb") as f:
                return f.read(3) == b"\xef\xbb\xbf"
        except OSError:
            return False

    def _decode_with_fallback_encodings(
        self, file_path: Path, original_error: Exception, kind: str
    ) -> Tuple[Optional[str], Optional[str], Optional[ValidationResult]]:
        """解析阶段的编码回退（审阅项 C-47 的另一半）。

        ``Bio.Phylo.parse`` / ``Bio.SeqIO.parse`` 拿到**路径**时按平台默认编码
        （UTF-8）读，非 UTF-8 输入会以原始 ``UnicodeDecodeError`` 冒到用户面前。
        这里按 :data:`TEXT_ENCODINGS` 真的重试一遍，把文本交回调用方解析。

        Returns:
            (text, encoding, None) 成功；(None, None, ValidationResult) 失败——
            失败文案给可执行的 ``iconv`` 建议，而不是"某个字节解不开"。
        """
        text, encoding = self._read_text_best_effort(file_path)
        if text is None:
            return (
                None,
                None,
                ValidationResult(
                    is_valid=False,
                    errors=[
                        f"{kind}文件不是 UTF-8，且按 {'/'.join(self.TEXT_ENCODINGS)} "
                        f"逐一重试后仍无法解码: {file_path}",
                        f"请用 iconv 转成 UTF-8 后重试，例如："
                        f"iconv -f GBK -t UTF-8 {file_path} > {file_path}.utf8"
                        f"（原始错误: {original_error}）",
                    ],
                ),
            )
        return text, encoding, None

    @staticmethod
    def _signature_matches(signature: Signature, data: bytes) -> bool:
        """比对/树格式判据的统一匹配方式（C-46 附带项）。

        判据表里的每一项可以是：

        * ``bytes`` —— 按**前缀**匹配（魔术字节/文件头）；
        * 已编译的 ``bytes`` 正则 —— 按**首行**匹配（PHYLIP 的
          "num_taxa seq_length" 头部行）。

        两种形态都由 :data:`TREE_FORMAT_SIGNATURES` /
        :data:`ALIGNMENT_FORMAT_SIGNATURES` 单一来源驱动，检测函数里不再各自
        硬编码一份判据。
        """
        if isinstance(signature, bytes):
            return data.startswith(signature)
        return signature.match(data) is not None

    def _phylip_header_matches(self, first_line: str) -> bool:
        """PHYLIP 头部行判据（取自 :data:`ALIGNMENT_FORMAT_SIGNATURES`）。

        格式检测与格式验证**共用这一条判据**，不再各写一遍正则——三处各写一份正是
        C-46 里"表是死数据、改表没有用"的成因。
        """
        data = first_line.encode("utf-8", errors="replace")
        return any(
            self._signature_matches(signature, data)
            for signature in self.ALIGNMENT_FORMAT_SIGNATURES.get(
                AlignmentFormat.PHYLIP, []
            )
        )

    @staticmethod
    def _looks_like_newick(text: str) -> bool:
        """Newick 的**词法**自检（审阅项 C-48：不依赖 BioPython）。

        判据保守：必须成对出现括号、至少两个末端标签（逗号或分号分隔）、
        不含 XML/Nexus 头。真正的结构校验留给 ``_try_parse_tree``。
        """
        stripped = text.strip()
        if not stripped or stripped[0] != "(":
            return False
        if stripped.count("(") != stripped.count(")"):
            return False
        if stripped.count("(") < 1:
            return False
        # 至少两个叶节点：单叶树 "A;" 不是可用输入
        if "," not in stripped and ";" not in stripped:
            return False
        upper = stripped[:10].upper()
        if upper.startswith("#NEXUS") or stripped.lstrip().startswith("<"):
            return False
        return True

    def _detect_tree_format(self, file_path: Path) -> Optional[TreeFormat]:
        """
        自动检测树文件格式

        Args:
            file_path: 树文件路径

        Returns:
            检测到的格式，无法检测返回 None
        """
        # 1. 根据扩展名判断（判据只从 TREE_FORMAT_BY_EXTENSION 读）
        detected = self.TREE_FORMAT_BY_EXTENSION.get(file_path.suffix.lower())
        if detected is not None:
            return detected

        # 2. 根据文件内容判断
        try:
            with open(file_path, "rb") as f:
                header = f.read(1024)  # 读取前 1KB
        except OSError as e:
            self.logger.warning(f"读取文件以检测树格式失败 ({file_path}): {e}")
            return None

        # 检查 Nexus / PhyloXML 格式（判据同样只从 TREE_FORMAT_SIGNATURES 读）
        for fmt in (TreeFormat.NEXUS, TreeFormat.PHYLOXML):
            for sig in self.TREE_FORMAT_SIGNATURES.get(fmt, []):
                if self._signature_matches(sig, header):
                    return fmt

        # Newick：先做与 BioPython 无关的词法自检（审阅项 C-48——环境里缺
        # biopython 时旧实现会得出"无法自动检测格式"，把矛头指向用户没问题的文件）
        text = header.decode("utf-8", errors="ignore")
        if self._looks_like_newick(text):
            return TreeFormat.NEWICK

        # 兜底：交给 BioPython 实解一次（不可用时显式说明是环境问题）
        try:
            from Bio import Phylo as BioPhylo
        except ImportError as e:
            self.logger.error(
                f"BioPython 不可用（{e}），无法进一步判定 {file_path} 的格式。"
                f"请执行 pip show biopython 确认已安装（本项目核心依赖）"
            )
            return None
        try:
            candidate = text.split(";")[0] + ";"
            if any(c in text for c in ["(", ")", ";", ":"]):
                BioPhylo.read(StringIO(candidate), "newick")
                return TreeFormat.NEWICK
        except Exception as e:
            self.logger.debug(f"内容兜底解析未通过，格式无法判定 ({file_path}): {e}")

        return None

    def _detect_alignment_format(self, file_path: Path) -> Optional[AlignmentFormat]:
        """
        自动检测序列比对文件格式

        Args:
            file_path: 序列比对文件路径

        Returns:
            检测到的格式，无法检测返回 None
        """
        # 1. 根据扩展名判断（判据只从 ALIGNMENT_FORMAT_BY_EXTENSION 读）
        detected = self.ALIGNMENT_FORMAT_BY_EXTENSION.get(file_path.suffix.lower())
        if detected is not None:
            return detected

        # 2. 根据文件内容判断
        try:
            with open(file_path, "rb") as f:
                header = f.read(4096)  # 读取前 4KB
        except OSError as e:
            # 审阅项 C-8：读不了文件是"检测失败"而不是"不是这些格式"，必须留痕
            self.logger.warning(f"读取文件以检测比对格式失败 ({file_path}): {e}")
            return None

        # 内容判据一律从 ALIGNMENT_FORMAT_SIGNATURES 读（C-46 附带项：旧的
        # 硬编码 b">" 与内联 PHYLIP 正则让这张表成了死数据）
        first_line = header.split(b"\n", 1)[0]
        for fmt in (
            AlignmentFormat.FASTA,
            AlignmentFormat.CLUSTAL,
            AlignmentFormat.STOCKHOLM,
            AlignmentFormat.PHYLIP,
        ):
            for signature in self.ALIGNMENT_FORMAT_SIGNATURES.get(fmt, []):
                # bytes 判据看整个文件头，正则判据只看首行（PHYLIP 的头部行）
                target = header if isinstance(signature, bytes) else first_line
                if self._signature_matches(signature, target):
                    return fmt

        return None

    def _validate_tree_format(
        self, file_path: Path, format: TreeFormat
    ) -> ValidationResult:
        """
        验证树文件格式

        审阅项 C-46：旧实现是 if/elif 分派且**没有 else**，于是
        ``TreeFormat.STOCKHOLDT`` 直接落到末尾的 ``return is_valid=True``——
        内容完全是垃圾文本的 ``junk.sth`` 也拿到"格式有效"。现在未覆盖的格式
        一律显式拒绝，而不是默认接受。

        审阅项 C-47：非 UTF-8 时**真的**按 utf-8-sig → utf-8 → cp1252 → latin-1
        依次重试，并在 warnings/info 里报告实际使用的编码；全部失败则给出可执行
        建议（iconv 转码），而不是把原始编解码报错丢给用户。

        Args:
            file_path: 树文件路径
            format: 期望的格式

        Returns:
            ValidationResult: 验证结果
        """
        warnings = []
        info = []

        content, encoding = self._read_text_best_effort(file_path, 4096)
        if content is None:
            return ValidationResult(
                is_valid=False,
                errors=[
                    f"无法以 {'/'.join(self.TEXT_ENCODINGS)} 任一编码读取 {file_path}",
                    "请用 iconv -f <源编码> -t UTF-8 <文件> 转码后重试；"
                    "常见来源是 Windows-1252 或带 BOM 的 UTF-8",
                ],
            )
        if encoding not in ("utf-8", "utf-8-sig"):
            warnings.append(
                f"文件 {file_path} 不是 UTF-8，已按 {encoding} 读取；"
                f"建议转成 UTF-8 以避免下游工具解读不一致"
            )
            # C-47：真正重试之后要把"用的到底是哪个编码"记进中性信息，
            # 便于用户/日志核对（旧的 warning 只承诺了一个不存在的回退行为）
            info.append(f"文本以 {encoding} 解码成功（{file_path}）")
        if encoding == "utf-8-sig":
            info.append(f"检测到 UTF-8 BOM（{file_path}）")

        if format == TreeFormat.NEWICK:
            # 检查 Newick 格式特征
            if not any(c in content for c in ["(", ")", ";"]):
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"文件不是有效的 Newick 格式: 缺少特征字符 '(', ')', ';': {file_path}"
                    ],
                )

        elif format == TreeFormat.NEXUS:
            if not content.lstrip().startswith("#NEXUS"):
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"文件不是有效的 Nexus 格式: 必须以 '#NEXUS' 开头: {file_path}"
                    ],
                )

        elif format == TreeFormat.PHYLOXML:
            if "<phyloxml" not in content.lower() and "<?xml" not in content.lower():
                return ValidationResult(
                    is_valid=False,
                    errors=[f"文件不是有效的 PhyloXML 格式: {file_path}"],
                )

        else:
            # STOCKHOLDT / AUTO 以及将来新增的枚举值：没做检查就不能声称通过（C-46）
            return ValidationResult(
                is_valid=False,
                format_detected=format.value,
                errors=[
                    f"未支持（本验证器对其不做任何实质检查）的树格式: "
                    f"{format.value} ({file_path})。树文件请用 Newick / Nexus / "
                    f"PhyloXML；Stockholm 是比对格式，不是树格式"
                ],
                warnings=warnings,
                info=info,
            )

        return ValidationResult(
            is_valid=True,
            format_detected=format.value,
            warnings=warnings,
            info=info,
        )

    def _validate_alignment_format(
        self, file_path: Path, format: AlignmentFormat
    ) -> ValidationResult:
        """
        验证序列比对文件格式

        审阅项 C-46：同 ``_validate_tree_format``——未覆盖的枚举值不再默认接受；
        ``PHYLIP_INTERLEAVED`` 此前**没有分支**（任意内容都被判"格式有效"），
        现在与 ``PHYLIP`` 共用同一套头部判据。

        Args:
            file_path: 序列比对文件路径
            format: 期望的格式

        Returns:
            ValidationResult: 验证结果
        """
        warnings = []
        info = []

        content, encoding = self._read_text_best_effort(file_path, 4096)
        if content is None:
            return ValidationResult(
                is_valid=False,
                errors=[
                    f"无法以 {'/'.join(self.TEXT_ENCODINGS)} 任一编码读取 {file_path}",
                    "请用 iconv -f <源编码> -t UTF-8 <文件> 转码后重试",
                ],
            )
        if encoding not in ("utf-8", "utf-8-sig"):
            warnings.append(
                f"文件 {file_path} 不是 UTF-8，已按 {encoding} 读取；建议转成 UTF-8"
            )
            # C-47：记录实际生效的编码（旧实现只留下"将尝试其他编码"的空头承诺）
            info.append(f"文本以 {encoding} 解码成功（{file_path}）")
        if encoding == "utf-8-sig":
            info.append(f"检测到 UTF-8 BOM（{file_path}）")

        if format == AlignmentFormat.FASTA:
            if not content.strip().startswith(">"):
                return ValidationResult(
                    is_valid=False,
                    errors=[f"文件不是有效的 FASTA 格式: 必须以 '>' 开头: {file_path}"],
                )

        elif format == AlignmentFormat.CLUSTAL:
            if not content.strip().upper().startswith("CLUSTAL"):
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"文件不是有效的 Clustal 格式: 必须以 'CLUSTAL' 开头: {file_path}"
                    ],
                )

        elif format == AlignmentFormat.STOCKHOLM:
            if not content.strip().upper().startswith("# STOCKHOLM"):
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"文件不是有效的 Stockholm 格式: 必须以 '# STOCKHOLM' 开头: {file_path}"
                    ],
                )

        elif format in (AlignmentFormat.PHYLIP, AlignmentFormat.PHYLIP_INTERLEAVED):
            # PHYLIP（简约式/交错式共用同一头部）：第一行是 "num_taxa seq_length"。
            # 判据取自 ALIGNMENT_FORMAT_SIGNATURES（C-46：此前这里是**第三份**内联
            # 重写的正则，表里的 PHYLIP 项因此形同死数据）
            stripped = content.strip()
            if not stripped:
                return ValidationResult(
                    is_valid=False, errors=[f"文件为空: {file_path}"]
                )
            first_line = stripped.split("\n", 1)[0].strip()
            if not self._phylip_header_matches(first_line):
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"文件不是有效的 PHYLIP 格式: 第一行应为 'num_taxa seq_length': {file_path}"
                    ],
                )

        else:
            return ValidationResult(
                is_valid=False,
                format_detected=format.value,
                errors=[
                    f"未支持（本验证器对其不做任何实质检查）的比对格式: "
                    f"{format.value} ({file_path})。请用 FASTA / PHYLIP / Clustal / "
                    f"Stockholm"
                ],
                warnings=warnings,
                info=info,
            )

        return ValidationResult(
            is_valid=True,
            format_detected=format.value,
            warnings=warnings,
            info=info,
        )

    def _try_parse_tree(self, file_path: Path, format: TreeFormat) -> ValidationResult:
        """
        尝试解析树文件

        Args:
            file_path: 树文件路径
            format: 树文件格式

        Returns:
            ValidationResult: 验证结果
        """
        warnings: List[str] = []
        info: List[str] = []

        try:
            from Bio import Phylo as BioPhylo
        except ImportError as e:
            # 审阅项 C-48：环境缺 BioPython 时旧实现把它和"文件坏了"压成同一条
            # "解析树文件失败"，报错把矛头指向用户本来没问题的文件。
            return ValidationResult(
                is_valid=False,
                errors=[
                    f"BioPython 不可用（{e}），无法完成解析校验: {file_path}",
                    "请先执行 pip show biopython 确认核心依赖已安装，"
                    "再重新运行（这不是树文件本身的问题）",
                ],
            )

        format_map = {
            TreeFormat.NEWICK: "newick",
            TreeFormat.NEXUS: "nexus",
            TreeFormat.PHYLOXML: "phyloxml",
        }

        bio_format = format_map.get(format)
        if bio_format is None:
            # C-46 同族：无法解析的格式不能声称"校验通过"
            return ValidationResult(
                is_valid=False,
                format_detected=format.value,
                errors=[
                    f"BioPython 无法解析 {format.value} 格式的树文件，"
                    f"因此无法确认其内容有效: {file_path}"
                ],
                warnings=warnings,
                info=info,
            )

        try:
            # 尝试解析所有树
            trees = list(BioPhylo.parse(str(file_path), bio_format))
        except UnicodeDecodeError as e:
            # 审阅项 C-47：旧实现在这里把原始编解码报错直接丢给用户
            # （"解析树文件失败: 'utf-8' codec can't decode byte 0xc9…"）。
            # 现在**真的**按 TEXT_ENCODINGS 重试一次并报告实际使用的编码；
            # 首选路径仍然是"把路径交给 Bio"（流式读，不把大文件整个读进内存），
            # 只有非 UTF-8 输入才付内存化的代价。
            text, encoding, failure = self._decode_with_fallback_encodings(
                file_path, original_error=e, kind="树"
            )
            if failure is not None:
                return failure
            try:
                trees = list(BioPhylo.parse(StringIO(text), bio_format))
            except Exception as retry_error:
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"解析树文件失败（已按 {encoding} 解码，但内容仍不是有效的"
                        f" {bio_format}）: {retry_error}"
                    ],
                )
            info.append(f"树文件以 {encoding} 解码后解析成功")
        except Exception as e:
            return ValidationResult(is_valid=False, errors=[f"解析树文件失败: {e}"])

        return self._summarize_parsed_tree(trees, file_path, warnings, info)

    def _summarize_parsed_tree(
        self,
        trees: List[Any],
        file_path: Path,
        warnings: List[str],
        info: List[str],
    ) -> ValidationResult:
        """把 Bio.Phylo 解析出的树列表汇总成验证结果（C-47 两条路径共用）。"""
        if not trees:
            return ValidationResult(
                is_valid=False, errors=[f"树文件中没有找到有效的树: {file_path}"]
            )

        tree_count = len(trees)

        # 检查第一棵树是否有叶节点
        tree = trees[0]
        terminals = list(tree.get_terminals())
        if len(terminals) < 2:
            return ValidationResult(
                is_valid=False,
                errors=[f"树文件中的树没有足够的叶节点 (至少需要 2 个): {file_path}"],
            )

        info.append(f"树文件包含 {len(terminals)} 个叶节点")

        # 检测多棵树
        if tree_count > 1:
            warnings.append(
                f"⚠ 树文件包含 {tree_count} 棵树！"
                f"请使用 --multi-tree-mode 参数指定处理方式：\n"
                f"    first  - 使用第一棵树\n"
                f"    last   - 使用最后一棵树\n"
                f"    random - 随机选择一棵树"
            )

        return ValidationResult(
            is_valid=True, warnings=warnings, info=info, tree_count=tree_count
        )

    def _try_parse_alignment(
        self, file_path: Path, format: AlignmentFormat
    ) -> ValidationResult:
        """
        尝试解析序列比对文件

        Args:
            file_path: 序列比对文件路径
            format: 序列比对文件格式

        Returns:
            ValidationResult: 验证结果
        """
        warnings: List[str] = []
        info: List[str] = []

        try:
            from Bio import SeqIO
        except ImportError as e:
            # 审阅项 C-48：与 _try_parse_tree 同理，缺依赖不该被说成"文件有问题"
            return ValidationResult(
                is_valid=False,
                errors=[
                    f"BioPython 不可用（{e}），无法完成解析校验: {file_path}",
                    "请先执行 pip show biopython 确认核心依赖已安装（这不是比对文件本身的问题）",
                ],
            )

        format_map = {
            AlignmentFormat.FASTA: "fasta",
            AlignmentFormat.PHYLIP: "phylip",
            AlignmentFormat.PHYLIP_INTERLEAVED: "phylip-sequential",
            AlignmentFormat.CLUSTAL: "clustal",
            AlignmentFormat.STOCKHOLM: "stockholm",
        }

        bio_format = format_map.get(format)
        if bio_format is None:
            return ValidationResult(
                is_valid=False,
                format_detected=format.value,
                errors=[
                    f"BioPython 无法解析 {format.value} 格式的比对文件，"
                    f"因此无法确认其内容有效: {file_path}"
                ],
                warnings=warnings,
                info=info,
            )

        try:
            # 尝试解析
            records = list(SeqIO.parse(str(file_path), bio_format))
        except UnicodeDecodeError as e:
            # 审阅项 C-47：与 _try_parse_tree 同一处理——真的按 TEXT_ENCODINGS
            # 重试，并把实际使用的编码报出来，而不是抛出原始编解码报错
            text, encoding, failure = self._decode_with_fallback_encodings(
                file_path, original_error=e, kind="序列比对"
            )
            if failure is not None:
                return failure
            try:
                records = list(SeqIO.parse(StringIO(text), bio_format))
            except Exception as retry_error:
                return ValidationResult(
                    is_valid=False,
                    errors=[
                        f"解析序列比对文件失败（已按 {encoding} 解码，但内容仍不是"
                        f"有效的 {bio_format}）: {retry_error}"
                    ],
                )
            info.append(f"序列比对文件以 {encoding} 解码后解析成功")
        except Exception as e:
            return ValidationResult(
                is_valid=False, errors=[f"解析序列比对文件失败: {e}"]
            )

        if not records:
            return ValidationResult(
                is_valid=False,
                errors=[f"序列比对文件中没有找到有效的序列: {file_path}"],
            )

        # 检查序列长度一致性（对于比对文件）
        if format != AlignmentFormat.FASTA:  # FASTA 不一定是比对后的
            lengths = set(len(record.seq) for record in records)
            if len(lengths) > 1:
                warnings.append(f"序列长度不一致，可能不是比对后的序列: {file_path}")

        info.append(f"序列比对文件包含 {len(records)} 条序列")

        return ValidationResult(is_valid=True, warnings=warnings, info=info)
