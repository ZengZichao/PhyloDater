"""
AlignmentMetadataExtractor - 多格式序列长度流式检测

支持 FASTA、PHYLIP、Stockholm、Clustal 格式的流式解析
"""

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from ..core.exceptions import AlignmentValidationError
from .logging import get_logger
from .safe_io import safe_writer


@dataclass
class AlignmentMetadata:
    """比对文件元数据"""

    format: str
    num_sequences: int
    sequence_length: int
    partition_lengths: Optional[List[int]] = None  # 多分区时各分区长度


class AlignmentMetadataExtractor:
    """
    比对元数据提取器

    使用流式解析避免将整个大文件加载到内存
    """

    MAX_INPUT_SIZE = 500 * 1024 * 1024

    def __init__(self) -> None:
        self.logger = get_logger()
        # CRLF 转换产生的临时副本路径；未转换时为 None。必须在 __init__ 里
        # 就把"可以为空"写进类型：旧代码只在 extract() 里随手赋 None，
        # 之后存 Path 就成了类型错。
        self._converted_path: Optional[Path] = None

    def extract(self, file_path: Path) -> AlignmentMetadata:
        """
        自动检测格式并提取元数据

        Args:
            file_path: 比对文件路径

        Returns:
            AlignmentMetadata: 比对元数据
        """
        if not file_path.exists():
            raise AlignmentValidationError(f"File not found: {file_path}")

        file_size = file_path.stat().st_size
        if file_size == 0:
            raise AlignmentValidationError(
                f"Empty file (0 bytes): {file_path}. "
                "Please provide a valid alignment file."
            )

        if file_size > self.MAX_INPUT_SIZE:
            raise AlignmentValidationError(
                f"File too large ({file_size / (1024*1024):.0f} MB), "
                f"exceeds limit ({self.MAX_INPUT_SIZE // (1024*1024)} MB). "
                "Use --low-memory mode or reduce data size."
            )

        # 每次 extract 都重置：上一份比对的临时副本不能泄到下一轮（实例可复用）。
        self._converted_path = None
        self._check_and_convert_line_endings(file_path)

        effective_path = self._converted_path or file_path

        try:
            fmt = self._detect_format(effective_path)

            if fmt == "fasta":
                return self._parse_fasta(effective_path)
            elif fmt == "phylip":
                return self._parse_phylip(effective_path)
            elif fmt == "stockholm":
                return self._parse_stockholm(effective_path)
            elif fmt == "clustal":
                return self._parse_clustal(effective_path)
            else:
                raise AlignmentValidationError(f"Unknown alignment format: {file_path}")
        finally:
            if self._converted_path and self._converted_path.exists():
                try:
                    self._converted_path.unlink()
                    self.logger.debug(
                        f"Cleaned up CRLF-converted temp file: {self._converted_path}"
                    )
                except Exception:
                    pass
                self._converted_path = None

    def _detect_format(self, file_path: Path) -> str:
        """检测比对文件格式"""
        ext = file_path.suffix.lower()

        # 根据扩展名初步判断
        if ext in [".fasta", ".fa", ".faa", ".fna", ".seq"]:
            return "fasta"
        elif ext in [".phylip", ".phy", ".ph"]:
            return "phylip"
        elif ext in [".stockholm", ".sto", ".stk"]:
            return "stockholm"
        elif ext in [".clustal", ".aln"]:
            return "clustal"

        # 根据首行内容判断
        with open(file_path, "r", encoding="utf-8", newline="") as f:
            first_line = f.readline()

        if first_line.startswith(">"):
            return "fasta"
        elif first_line.startswith("# STOCKHOLM"):
            return "stockholm"
        elif first_line.startswith("CLUSTAL") or "CLUSTAL" in first_line:
            return "clustal"
        elif re.match(r"^\s*\d+\s+\d+", first_line):
            return "phylip"

        return "unknown"

    def _parse_fasta(self, file_path: Path) -> AlignmentMetadata:
        """流式解析 FASTA 格式"""
        num_seqs = 0
        lengths = []
        current_length = 0

        with open(file_path, "r", encoding="utf-8", newline="") as f:
            for line in f:
                line = line.strip()

                if line.startswith(">"):
                    # 新序列开始
                    if num_seqs > 0:
                        lengths.append(current_length)
                        # 前三条序列长度不一致则提前报错
                        if len(lengths) >= 3 and len(set(lengths)) > 1:
                            raise AlignmentValidationError(
                                f"Sequence length mismatch in FASTA: {lengths}"
                            )
                    num_seqs += 1
                    current_length = 0
                elif line:
                    # 序列内容
                    current_length += len(line)

            # 最后一条序列
            if num_seqs > 0:
                lengths.append(current_length)

        if num_seqs == 0:
            raise AlignmentValidationError("No sequences found in FASTA file")

        # 验证所有序列长度一致
        unique_lengths = set(lengths)
        if len(unique_lengths) > 1:
            raise AlignmentValidationError(
                f"Sequence length mismatch: {sorted(unique_lengths)}"
            )

        return AlignmentMetadata(
            format="fasta", num_sequences=num_seqs, sequence_length=lengths[0]
        )

    def _parse_phylip(self, file_path: Path) -> AlignmentMetadata:
        """解析 PHYLIP 格式"""
        with open(file_path, "r", encoding="utf-8", newline="") as f:
            # 读取首行
            first_line = f.readline().strip()
            match = re.match(r"^\s*(\d+)\s+(\d+)", first_line)

            if not match:
                raise AlignmentValidationError("Invalid PHYLIP header")

            num_seqs = int(match.group(1))
            seq_len = int(match.group(2))

            # 检查是否为多分区格式 (G 选项)
            if "G" in first_line.upper():
                # 读取分区信息
                partition_line = f.readline().strip()
                parts = partition_line.split()
                if len(parts) >= 2:
                    num_partitions = int(parts[0])
                    partition_lengths = [int(x) for x in parts[1 : num_partitions + 1]]

                    return AlignmentMetadata(
                        format="phylip",
                        num_sequences=num_seqs,
                        sequence_length=sum(partition_lengths),
                        partition_lengths=partition_lengths,
                    )

        return AlignmentMetadata(
            format="phylip", num_sequences=num_seqs, sequence_length=seq_len
        )

    def _parse_stockholm(self, file_path: Path) -> AlignmentMetadata:
        """流式解析 Stockholm 格式"""
        num_seqs = 0
        seq_lengths = {}

        with open(file_path, "r", encoding="utf-8", newline="") as f:
            for line in f:
                line = line.strip()

                # 跳过注释行
                if line.startswith("#") or line.startswith("//"):
                    continue

                # 序列行
                parts = line.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    seq = parts[1]

                    if seq_id not in seq_lengths:
                        num_seqs += 1
                        seq_lengths[seq_id] = 0

                    seq_lengths[seq_id] += len(seq)

        if num_seqs == 0:
            raise AlignmentValidationError("No sequences found in Stockholm file")

        # 验证长度一致
        lengths = list(seq_lengths.values())
        unique_lengths = set(lengths)
        if len(unique_lengths) > 1:
            raise AlignmentValidationError(
                f"Sequence length mismatch: {sorted(unique_lengths)}"
            )

        return AlignmentMetadata(
            format="stockholm", num_sequences=num_seqs, sequence_length=lengths[0]
        )

    def _parse_clustal(self, file_path: Path) -> AlignmentMetadata:
        """流式解析 Clustal 格式"""
        num_seqs = 0
        seq_lengths = {}

        with open(file_path, "r", encoding="utf-8", newline="") as f:
            # 跳过首行（CLUSTAL 标识）
            next(f)

            for line in f:
                line = line.strip()

                if not line or line.startswith("//"):
                    continue

                # 序列行
                parts = line.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    seq = parts[1]

                    if seq_id not in seq_lengths:
                        num_seqs += 1
                        seq_lengths[seq_id] = 0

                    seq_lengths[seq_id] += len(seq)

        if num_seqs == 0:
            raise AlignmentValidationError("No sequences found in Clustal file")

        # 验证长度一致
        lengths = list(seq_lengths.values())
        unique_lengths = set(lengths)
        if len(unique_lengths) > 1:
            raise AlignmentValidationError(
                f"Sequence length mismatch: {sorted(unique_lengths)}"
            )

        return AlignmentMetadata(
            format="clustal", num_sequences=num_seqs, sequence_length=lengths[0]
        )

    def validate_sequence_type(self, file_path: Path) -> str:
        """
        自动检测序列类型（DNA 或蛋白质）

        使用计数策略：统计 DNA-only 字符和 protein-only 字符的比例。
        IUPAC 模糊码（B, D, H, V, Z 等）同时出现在 DNA 和蛋白质中，
        因此不作为判据。仅当 protein-only 字符显著超过 DNA-only 字符时
        才判定为蛋白质。

        Returns:
            'dna' 或 'protein'
        """
        # DNA 独有字符（不出现在标准蛋白质字母表中）
        dna_only = set("ATCGU")
        # 蛋白质独有字符（不出现在 IUPAC DNA 模糊码中）
        protein_only = set("EFILPQSVWY")
        # 可疑字符：B, D, H, R, N, Z 等在 DNA 和蛋白质中均可能出现
        # 不计入任一类

        dna_count = 0
        protein_count = 0

        with open(file_path, "r", encoding="utf-8", newline="") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">") or line.startswith("#"):
                    continue
                for char in line.upper():
                    if char in dna_only:
                        dna_count += 1
                    elif char in protein_only:
                        protein_count += 1
                    if dna_count + protein_count >= 1000:
                        break
                if dna_count + protein_count >= 1000:
                    break

        # 仅当 protein-only 字符超过总量的 5% 且显著多于 DNA-only 时判定为蛋白质
        total = dna_count + protein_count
        if total == 0:
            return "dna"
        return (
            "protein"
            if (protein_count > total * 0.05 and protein_count > dna_count)
            else "dna"
        )

    def _check_and_convert_line_endings(self, file_path: Path) -> None:
        """
        检测 Windows 换行符 (CRLF) 并创建转换后的副本

        不修改用户原始文件，而是在同目录下创建 .lf 转换副本。
        返回转换后文件的路径（如果无需转换则返回原路径）。
        """
        # 先检查文件是否包含 CRLF（只读前 8KB）
        with open(file_path, "rb") as f:
            sample = f.read(8192)

        if b"\r\n" not in sample:
            return  # 快速路径：无 CRLF

        self.logger.warning(
            f"Detected Windows line endings (CRLF) in {file_path.name}. "
            "Creating a converted copy with Unix line endings (LF)."
        )

        # 流式转换：读取整个文件并替换
        with open(file_path, "rb") as f:
            raw_content = f.read()

        unix_content = raw_content.replace(b"\r\n", b"\n")

        # 写入系统临时目录，避免在用户数据目录留下文件
        import tempfile

        try:
            fd, tmp_name = tempfile.mkstemp(suffix=".lf", prefix="phylodater_")
            os.close(fd)
            converted_path = Path(tmp_name)
            with safe_writer(converted_path, "wb") as f:
                f.write(unix_content)
            os.chmod(converted_path, 0o600)
        except OSError as e:
            self.logger.warning(
                f"无法创建 LF 转换副本（{e}），将继续使用原文件，"
                "若解析失败请手动转换为 Unix 换行。"
            )
            return

        self.logger.info(f"LF-converted copy saved to: {converted_path.name}")

        # 用转换后的文件路径替换原始路径
        # 注意：调用方应使用转换后的路径
        self._converted_path = converted_path

    def validate_alignment_file(self, file_path: Path) -> bool:
        """
        验证比对文件的完整性

        检查项目：
        1. 文件存在
        2. 非空
        3. 有效的行结尾

        Returns:
            True 如果验证通过
        """
        if not file_path.exists():
            raise AlignmentValidationError(f"File not found: {file_path}")

        file_size = file_path.stat().st_size
        if file_size == 0:
            raise AlignmentValidationError(
                f"Empty file (0 bytes): {file_path}. "
                "Please provide a valid alignment file."
            )

        with open(file_path, "rb") as f:
            raw = f.read(4096)

        if b"\r\n" in raw:
            self.logger.warning(
                f"File {file_path.name} contains Windows line endings (CRLF). "
                "This may cause parsing errors. Consider converting to Unix line endings."
            )

        return True
