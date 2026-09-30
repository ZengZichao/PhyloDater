"""
PhyloDater 适配器层

包含各定年软件的适配器实现
"""

from .lsd2_method import LSD2Method
from .mcmctree_method import MCMCTreeMethod
from .mdcat_method import MDCatMethod

# 导入适配器以自动注册
from .pathd8_method import PATHd8Method
from .r8s_pyr8s_method import R8sPyr8sMethod
from .treepl_method import TreePLMethod
from .wlogdate_method import WLogDateMethod

__all__ = [
    "PATHd8Method",
    "R8sPyr8sMethod",
    "TreePLMethod",
    "LSD2Method",
    "WLogDateMethod",
    "MCMCTreeMethod",
    "MDCatMethod",
]
