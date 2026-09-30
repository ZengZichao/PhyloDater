"""集中式落盘写入入口：路径穿越纵深防御。

所有对用户可见路径的写文件操作都经由 :func:`safe_writer`：

1. 拒绝任何含 ``..`` 段的路径——防止经由配置值拼出逃逸工作目录的路径；
2. ``resolve()`` 规范化符号链接与相对段后再落盘。

底层用 ``os.open`` / ``os.fdopen`` 实现，打开语义（``O_WRONLY|O_CREAT|O_TRUNC``、
``0666 & ~umask``、文本模式的 ``encoding``/``newline`` 等参数透传）与内置
``open`` 保持一致；唯一的可见差异是含 ``..`` 段的路径会抛
``ValueError`` 而不再静默写出目标目录之外。
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

__all__ = ["safe_writer"]


def safe_writer(path: Path | str, mode: str = "w", **kwargs: Any) -> Any:
    """以 ``mode`` 打开 ``path`` 供写入；含 ``..`` 段的路径直接拒绝。"""
    target = Path(path)
    if ".." in target.parts:
        raise ValueError(f"path must not contain '..' segments: {target}")
    resolved = target.resolve()
    fd = os.open(resolved, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o666)
    try:
        return os.fdopen(fd, mode, **kwargs)
    except Exception:
        os.close(fd)
        raise
