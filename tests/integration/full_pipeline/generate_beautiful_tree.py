#!/usr/bin/env python3
"""
生成更美观的时间树图片
"""

import sys
from pathlib import Path

# 项目根目录
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ete3 import Tree
from matplotlib.collections import LineCollection

from phylodater.viz.tree_plot import compute_node_coordinates, get_tree_segments

# 创建测试树 - 使用更复杂的树
newick = """(
    (
        (
            (human:0.007, chimpanzee:0.007):0.006,
            gorilla:0.013
        ):0.009,
        (
            orangutan:0.018,
            (gibbon:0.022, siamang:0.022):0.004
        ):0.015
    ):0.011,
    (
        (
            (macaque:0.025, baboon:0.025):0.012,
            (vervet:0.032, colobus:0.032):0.008
        ):0.007,
        (
            (lemur:0.045, loris:0.045):0.010,
            tarsier:0.050
        ):0.015
    ):0.010
):0.010;"""

# 创建ete3树
ete_tree = Tree(newick, format=1)

# 计算坐标 - 使用55 Ma作为根年龄
root_age = 55.0
coords = compute_node_coordinates(ete_tree, root_age=root_age)
segments, properties = get_tree_segments(ete_tree, coords, root_age)

# 创建图片
fig, (ax_tree, ax_geo) = plt.subplots(
    2, 1, figsize=(16, 10), gridspec_kw={"height_ratios": [0.85, 0.15], "hspace": 0.05}
)

# ==================== 绘制树 ====================

# 定义颜色方案
BRANCH_COLOR = "#2C3E50"
LEAF_COLOR = "#E74C3C"
INTERNAL_COLOR = "#3498DB"
ROOT_COLOR = "#27AE60"

# 绘制分支线段
lc = LineCollection(segments, colors=BRANCH_COLOR, linewidths=2.0, zorder=1, alpha=0.8)
ax_tree.add_collection(lc)

# 绘制节点
leaf_x, leaf_y = [], []
internal_x, internal_y = [], []
root_x, root_y = [], []

for node, coord in coords.items():
    if coord.is_root:
        root_x.append(coord.x)
        root_y.append(coord.y)
    elif coord.is_leaf:
        leaf_x.append(coord.x)
        leaf_y.append(coord.y)
    else:
        internal_x.append(coord.x)
        internal_y.append(coord.y)

# 绘制叶节点
if leaf_x:
    ax_tree.scatter(
        leaf_x,
        leaf_y,
        s=80,
        c=LEAF_COLOR,
        zorder=3,
        edgecolors="white",
        linewidths=1.5,
        label="Leaf nodes",
    )

# 绘制内部节点
if internal_x:
    ax_tree.scatter(
        internal_x,
        internal_y,
        s=40,
        c=INTERNAL_COLOR,
        zorder=3,
        edgecolors="white",
        linewidths=1.0,
        alpha=0.7,
        label="Internal nodes",
    )

# 绘制根节点
if root_x:
    ax_tree.scatter(
        root_x,
        root_y,
        s=100,
        c=ROOT_COLOR,
        zorder=3,
        edgecolors="white",
        linewidths=2.0,
        marker="*",
        label="Root",
    )

# 添加叶节点标签
for node, coord in coords.items():
    if coord.is_leaf and coord.name:
        ax_tree.annotate(
            coord.name,
            xy=(coord.x, coord.y),
            xytext=(5, 0),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
            color="#2C3E50",
            va="center",
            fontstyle="italic",
        )

# 添加年龄标签到内部节点
for node, coord in coords.items():
    if not coord.is_leaf and not coord.is_root and coord.x > 5:
        ax_tree.annotate(
            f"{coord.x:.1f}",
            xy=(coord.x, coord.y),
            xytext=(0, 12),
            textcoords="offset points",
            fontsize=7,
            color="#7F8C8D",
            ha="center",
            style="italic",
        )

# 设置树的坐标轴
ax_tree.set_xlim(-2, root_age + 5)
ax_tree.set_ylim(0, len(leaf_x) + 1)
ax_tree.invert_xaxis()  # 年龄从右到左

# 添加网格
ax_tree.grid(axis="x", alpha=0.3, linestyle="--", color="#BDC3C7")
ax_tree.set_axisbelow(True)

# 设置标签
ax_tree.set_ylabel("Taxa", fontsize=12, fontweight="bold")
ax_tree.set_xlabel("")  # 不显示x轴标签，因为有地质年代

# 添加图例
ax_tree.legend(loc="upper left", fontsize=9, framealpha=0.9, edgecolor="#BDC3C7")

# 添加标题
ax_tree.set_title(
    "Phylogenetic Time Tree (MCMCTree)", fontsize=16, fontweight="bold", pad=15
)

# ==================== 绘制地质年代 ====================

# 定义地质年代数据 (简化版)
geo_periods = [
    # (名称, 开始年龄, 结束年龄, 颜色)
    ("Holocene", 0, 0.012, "#FFFF00"),
    ("Pleistocene", 0.012, 2.58, "#FFE4B5"),
    ("Pliocene", 2.58, 5.33, "#FFA500"),
    ("Miocene", 5.33, 23.03, "#FFD700"),
    ("Oligocene", 23.03, 33.9, "#90EE90"),
    ("Eocene", 33.9, 56.0, "#00CED1"),
]

# 绘制地质年代条
for name, start, end, color in geo_periods:
    if end <= root_age:
        ax_geo.axhspan(
            0, 1, xmin=start / root_age, xmax=end / root_age, color=color, alpha=0.7
        )
        # 添加标签
        mid = (start + end) / 2
        if mid < root_age:
            ax_geo.text(
                mid,
                0.5,
                name,
                ha="center",
                va="center",
                fontsize=8,
                fontweight="bold",
                rotation=0,
            )

# 设置地质年代坐标轴
ax_geo.set_xlim(ax_tree.get_xlim())
ax_geo.set_ylim(0, 1)
ax_geo.invert_xaxis()
ax_geo.set_yticks([])
ax_geo.set_xlabel("Age (Ma)", fontsize=12, fontweight="bold")

# 添加年代刻度
ax_geo.xaxis.set_major_locator(plt.MultipleLocator(10))
ax_geo.xaxis.set_minor_locator(plt.MultipleLocator(5))

# ==================== 美化 ====================

# 设置整体样式
for ax in [ax_tree, ax_geo]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(width=0.5, labelsize=10)

# 添加水印
fig.text(
    0.98,
    0.02,
    "Generated by PhyloDater",
    fontsize=8,
    color="#BDC3C7",
    ha="right",
    va="bottom",
    style="italic",
    alpha=0.7,
)

# 保存图片到项目内 output 目录
OUTPUT_DIR = PROJECT_ROOT / "output" / "integration_pipeline"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
output_path = OUTPUT_DIR / "phylodater_timetree_beautiful.png"
plt.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
    pad_inches=0.2,
)

print(f"图片已保存到: {output_path}")
print(f'文件大小: {__import__("os").path.getsize(output_path)} bytes')
