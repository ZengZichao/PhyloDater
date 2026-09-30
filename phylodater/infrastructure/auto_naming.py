"""
Auto Naming - 自动目录命名功能

根据输入文件自动生成输出目录名称
"""

import re
from datetime import datetime
from pathlib import Path
from typing import Optional

from .logging import get_logger


class AutoNamingManager:
    """
    自动命名管理器

    根据输入文件名自动生成输出目录名称
    """

    def __init__(self) -> None:
        self.logger = get_logger()

    def generate_output_dir(
        self,
        input_file: Path,
        method: str,
        suffix: Optional[str] = None,
        timestamp: bool = True,
    ) -> Path:
        """
        生成输出目录名称

        格式: {input_stem}_{method}_{suffix}_{timestamp}

        Args:
            input_file: 输入文件路径
            method: 方法名称 (如 lsd2, mcmctree)
            suffix: 可选后缀
            timestamp: 是否添加时间戳
        """
        # 获取输入文件的基本名称
        input_stem = input_file.stem

        # 清理名称（移除特殊字符）
        input_stem = self._sanitize_name(input_stem)

        # 构建目录名
        parts = [input_stem, method]

        if suffix:
            suffix = self._sanitize_name(suffix)
            parts.append(suffix)

        if timestamp:
            timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
            parts.append(timestamp_str)

        dir_name = "_".join(parts)

        # 确保目录名合法
        dir_name = self._ensure_valid_dirname(dir_name)

        output_dir = input_file.parent / dir_name

        self.logger.info(f"Auto-generated output directory: {output_dir}")

        return output_dir

    def generate_run_dir(self, base_dir: Path, run_id: int, method: str) -> Path:
        """
        生成运行子目录名称

        Args:
            base_dir: 基础目录
            run_id: 运行ID
            method: 方法名称
        """
        dir_name = f"{method}_run{run_id}"
        return base_dir / dir_name

    def _sanitize_name(self, name: str) -> str:
        """清理名称，替换特殊字符"""
        # 替换特殊字符为下划线
        sanitized = re.sub(r"[^\w\-]", "_", name)
        # 移除连续的下划线
        sanitized = re.sub(r"_+", "_", sanitized)
        # 移除首尾下划线
        sanitized = sanitized.strip("_")
        return sanitized

    def _ensure_valid_dirname(self, name: str) -> str:
        """确保目录名合法"""
        # 限制长度
        max_length = 100
        if len(name) > max_length:
            name = name[:max_length]

        # 确保不以点开头（隐藏文件）
        if name.startswith("."):
            name = "run_" + name

        # 确保不为空
        if not name:
            name = "output"

        return name

    def find_unique_dir(self, base_dir: Path) -> Path:
        """
        如果目录已存在，添加数字后缀使其唯一

        Args:
            base_dir: 基础目录路径
        """
        if not base_dir.exists():
            return base_dir

        counter = 1
        while True:
            new_dir = Path(f"{base_dir}_{counter}")
            if not new_dir.exists():
                self.logger.info(
                    f"Directory {base_dir} exists, using {new_dir} instead"
                )
                return new_dir
            counter += 1

            # 防止无限循环
            if counter > 1000:
                raise RuntimeError(f"Cannot find unique directory name for {base_dir}")


class LSD2AutoNaming(AutoNamingManager):
    """LSD2 专用自动命名"""

    def generate_from_alignment_tree(
        self, alignment_file: Path, tree_file: Path, output_dir: Optional[Path] = None
    ) -> Path:
        """
        根据比对和树文件生成输出目录

        参考 run_lsd2.py 的命名逻辑
        """
        if output_dir:
            return output_dir

        # 使用树文件名作为基础
        tree_stem = tree_file.stem

        # 移除常见后缀
        tree_stem = re.sub(r"_(tree|tre|nwk|newick)$", "", tree_stem, flags=re.I)

        return self.generate_output_dir(
            tree_file, method="lsd2", suffix="dating", timestamp=True
        )


class MCMCTreeAutoNaming(AutoNamingManager):
    """MCMCTree 专用自动命名"""

    def generate_pipeline_dir(
        self, seq_file: Path, tree_file: Path, output_base: Optional[str] = None
    ) -> Path:
        """
        生成 MCMCTree 流程目录

        参考 mcmctree-toolkit 的命名逻辑
        """
        if output_base:
            base_path = Path(output_base)
        else:
            # 使用序列文件名
            base_path = seq_file.parent / f"{seq_file.stem}_mcmctree"

        return self.find_unique_dir(base_path)
