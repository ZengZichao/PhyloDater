"""PhyloDater 版本信息。

单一事实来源：版本号同时写在此处与 ``pyproject.toml`` 的 ``[project] version``，
二者必须一致。此处刻意不使用 git 派生的动态版本（无 ``.devN``/提交哈希等迭代串）。
"""

from __future__ import annotations

__all__ = [
    "__version__",
    "__version_tuple__",
    "version",
    "version_tuple",
]

__version__ = version = "0.1.0"
__version_tuple__ = version_tuple = (0, 1, 0)
