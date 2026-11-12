# -*- coding: utf-8 -*-
"""绘制论文图5：包装器正确性对照（a）与语义退化信息损失评估（b）。

数据来源（PhyloDater 项目代码仓库，只读）：
- tier1_manual_results.json   手工配置原生工具运行 vs PhyloDater 包装器
- degradation_results.json    MCMCTree gamma 先验 → 95% CI 均匀分布退化的节点年龄偏差

输出：与脚本同目录下的 图5-包装正确性与语义退化评估.pdf / .png

本轮修订要点（纯视觉精修：数据编码、面板划分与图注承诺全部不变）：
- 版面：改为手工 add_axes 定位 + bbox_inches="tight"。两栏共用同一 y 区间与高度，
  因此标题行与图例行天然对齐；原先 (b) 的横轴标签与标题被画布裁掉的问题一并消除。
- (a) 纵轴由 ±0.05 Ma 收紧到 ±0.015 Ma：容差带从只占画面 1/5 变为占据主体，
  大面积空白消失，"19 个点全部落在带内"这一结论一眼可读。
- (a) 容差带补 ±0.01 Ma 边界细线（与纵轴刻度严格对齐），带体内加 0/5/10/15 纵向
  刻度线以便按原生年龄读数；一致线改为更细的虚线，不再压过标记。
- (a) 标记加白色描边、beeswarm 最小间距由 0.40 提到 0.55 Ma（略大于标记直径），
  密集处不再糊成一团。
- 图例：两栏统一改为"轴上方单行居中"的页眉式图例。(a) 不再用竖排图例悬在右上
  空白里，(b) 不再与 HumanChimp 行的标记和数值标签叠字。
- (b) 横轴左端由 −6.5 扩到 −8.0、行名与绘图区留白加大，数值标签不再撞到行名；
  加浅灰隔行底纹与纵向网格，把远离绘图区的行名与短棒在视觉上连起来；数值标签改用
  对应系列色，标签与系列直接绑定；棒端改圆头。
- 体例：面板标题改为与图 1-4 一致的深板岩色粗体、左对齐于轴左边、"(a)" 后单空格；
  轴线、刻度、坐标轴标签与正文文字统一为同一套板岩灰阶，负号一律 U+2212。
"""

import json
import os
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# 运行时自动切换工作目录到脚本所在路径
SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

# 动态获取脚本名（不含路径与扩展名），禁止硬编码
SCRIPT_NAME = Path(sys.argv[0]).resolve().stem if sys.argv[0] else Path(__file__).stem


def _find_bench():
    """定位 examples/benchmark/。
    仓库内副本（scripts/make_manuscript_figures.py）以仓库根为基准解析：
    本脚本位于 <repo>/scripts/，``examples/benchmark`` 是它的同级目录。
    """
    cands = []
    env = os.environ.get("PHYLODATER_REPO")
    if env:
        cands.append(Path(env) / "examples" / "benchmark")
    here = SCRIPT_DIR
    repo_root = here.parent
    cands += [
        repo_root / "examples" / "benchmark",
        here / "examples" / "benchmark",
        here.parent / "PhyloDater-项目代码" / "examples" / "benchmark",
    ]
    for c in cands:
        if (c / "tier1_manual_results.json").is_file():
            return c
    raise SystemExit(
        "examples/benchmark with tier1_manual_results.json and "
        "degradation_results.json not found; set PHYLODATER_REPO to the "
        "PhyloDater repository root."
    )


BENCH = _find_bench()
TIER1_JSON = BENCH / "tier1_manual_results.json"
DEG_JSON = BENCH / "degradation_results.json"

MINUS = "\u2212"  # U+2212 typographic minus, aligned with axis ticks

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Helvetica",
        "mathtext.it": "Helvetica:italic",
        "mathtext.bf": "Helvetica:bold",
        "mathtext.default": "regular",
        "pdf.fonttype": 42,  # TrueType instead of Type 3 (journal requirement)
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "axes.linewidth": 0.6,
        "axes.axisbelow": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.color": "#5B7078",
        "ytick.color": "#5B7078",
        "axes.edgecolor": "#5B7078",
        "axes.labelcolor": "#2C4149",
        "text.color": "#2C4149",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

# 由 Okabe-Ito 降饱和、统一明度后得到的一组配色：与图 1-4 的石油蓝绿/鼠尾草绿/
# 砂土黄/砖玫红体系同族，印刷时不与中性底色串色
METHOD_STYLE = {
    "pathd8": dict(color="#2C6E86", marker="o", label="PATHd8"),
    "lsd2": dict(color="#BE6B2C", marker="s", label="LSD2"),
    "treepl": dict(color="#3F8A6C", marker="^", label="treePL"),
    "wlogdate": dict(color="#9C5F87", marker="D", label="wLogDate"),
}
SEED_STYLE = {
    "42": dict(color="#2C6E86", marker="o", label="seed 42"),
    "43": dict(color="#BE6B2C", marker="s", label="seed 43"),
}
NODE_ROWS = [
    ("HumanChimp", "HumanChimp"),
    ("Root", "GreatApes = Root"),
    ("internal_node_1", "internal node 1"),
    ("internal_node_2", "internal node 2"),
    ("internal_node_3", "internal node 3"),
]

# 中性色：与图 1-4 的浅灰底、板岩描文同族
INK = "#2C4149"  # 标题与正文
SLATE = "#5B7078"  # 轴线与参考线
GRID = "#DCE3E3"  # 纵向网格
ROW_TINT = "#F3F6F6"  # 隔行底纹
BAND_FILL = "#E7EEED"  # ±0.01 Ma 容差带
BAND_EDGE = "#B3C3C3"  # 容差带边界
GRID_A = "#CFDADA"  # (a) 纵向网格，压在带体上

TOLERANCE_MA = 0.01  # Ma —— 图注赖以立论的容差带
MARKER_MIN_GAP_MA = 0.62  # 标记直径 + 白描边，相邻标记之间留出可见缝隙
MARKER_SIZE = 26


def beeswarm(xs, min_gap):
    """确定性一维 beeswarm：仅在重叠处按最小间距推开，其余保持原横坐标。"""
    xs = np.asarray(xs, dtype=float)
    order = np.argsort(xs, kind="stable")
    v = xs[order].copy()
    n = v.size
    if n < 2:
        return xs.copy()
    # 1) 左到右强制满足最小间距
    for i in range(1, n):
        if v[i] - v[i - 1] < min_gap:
            v[i] = v[i - 1] + min_gap
    # 2) 双向松弛：在邻居约束下把每个点拉回其真实横坐标
    for _ in range(400):
        moved = False
        for i in range(n - 1, -1, -1):
            hi = v[i + 1] - min_gap if i < n - 1 else np.inf
            want = min(xs[order][i], hi)
            if want < v[i] - 1e-12:
                v[i], moved = want, True
        for i in range(n):
            lo = v[i - 1] + min_gap if i > 0 else -np.inf
            want = max(xs[order][i], lo)
            if want > v[i] + 1e-12:
                v[i], moved = want, True
        if not moved:
            break
    # 3) 对每个紧贴成串的簇做整体居中，消除单向漂移
    i = 0
    while i < n:
        j = i
        while j + 1 < n and v[j + 1] - v[j] <= min_gap + 1e-9:
            j += 1
        block = v[i : j + 1]
        shift = xs[order][i : j + 1].mean() - block.mean()
        left_ok = (v[i - 1] + min_gap) - block[0] if i > 0 else -np.inf
        right_ok = (v[j + 1] - min_gap) - block[-1] if j < n - 1 else np.inf
        v[i : j + 1] = block + float(np.clip(shift, left_ok, right_ok))
        i = j + 1
    out = np.empty_like(xs)
    out[order] = v
    return out


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _header_legend(ax, ncol, columnspacing, handlelength, handletextpad):
    """图例统一放在轴上方居中，与标题共组成两行页眉，避免压住数据。"""
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.015),
        ncol=ncol,
        frameon=False,
        fontsize=8,
        borderaxespad=0,
        columnspacing=columnspacing,
        handlelength=handlelength,
        handletextpad=handletextpad,
    )


def panel_a(ax, tier1):
    """(a) 残差图：wrapper − native vs native，容差带 ±0.01 Ma。"""
    style, manual, residual = [], [], []
    for method, info in tier1["methods"].items():
        s = METHOD_STYLE[method]
        for row in info["comparison"].values():
            if row["wrapper"] is None or row["manual"] is None:
                continue
            style.append(s)
            manual.append(row["manual"])
            residual.append(row["wrapper"] - row["manual"])
    manual = np.asarray(manual, dtype=float)
    residual = np.asarray(residual, dtype=float)
    plotted_x = beeswarm(manual, MARKER_MIN_GAP_MA)

    x_lo, x_hi = -0.9, 15.9
    y_lim = 0.015

    # 容差带：带体 + 与刻度严格对齐的边界细线
    ax.fill_between(
        [x_lo, x_hi],
        -TOLERANCE_MA,
        TOLERANCE_MA,
        color=BAND_FILL,
        zorder=1,
        linewidth=0,
    )
    # 纵向刻度线压过带体，便于按原生年龄读数
    for gx in (0, 5, 10, 15):
        ax.vlines(gx, -TOLERANCE_MA, TOLERANCE_MA, color=GRID_A, lw=0.6, zorder=1.5)
    ax.hlines(
        [TOLERANCE_MA, -TOLERANCE_MA], x_lo, x_hi, color=BAND_EDGE, lw=0.55, zorder=2
    )
    # 完全一致线
    ax.axhline(0.0, color=SLATE, lw=0.85, ls=(0, (3.4, 2.6)), zorder=3)

    for s in METHOD_STYLE.values():
        ax.scatter(
            [],
            [],
            s=MARKER_SIZE,
            color=s["color"],
            marker=s["marker"],
            edgecolors="white",
            linewidths=0.5,
            label=s["label"],
        )
    for s, x, y in zip(style, plotted_x, residual):
        ax.scatter(
            x,
            y,
            s=MARKER_SIZE,
            color=s["color"],
            marker=s["marker"],
            edgecolors="white",
            linewidths=0.5,
            zorder=5,
        )

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-y_lim, y_lim)
    ax.set_xticks([0, 5, 10, 15])
    ax.set_yticks([-TOLERANCE_MA, 0.0, TOLERANCE_MA], [MINUS + "0.01", "0", "+0.01"])
    ax.set_xlabel("Manually configured native run (Ma)", labelpad=4)
    ax.set_ylabel("Wrapper \u2212 native (Ma)", labelpad=4)
    ax.tick_params(axis="y", which="major", pad=3)
    ax.text(
        0.99,
        0.94,
        f"n = {len(manual)} node\u2013method pairs",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color=SLATE,
    )
    _header_legend(
        ax, ncol=4, columnspacing=1.05, handlelength=0.95, handletextpad=0.35
    )


def _signed_pct(v):
    """带符号百分比，负号用 U+2212。"""
    s = f"{abs(v):.1f}%"
    return (MINUS + s) if v < 0 else ("+" + s)


def panel_b(ax, deg):
    """(b) 带零线棒图：节点年龄相对 gamma 先验的符号偏差（%）。"""
    seeds = ["42", "43"]
    y = np.arange(len(NODE_ROWS))
    x_lo, x_hi = -8.0, 11.2

    # 隔行底纹：把远离绘图区的行名与短棒在视觉上连起来
    for i in range(len(NODE_ROWS)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=ROW_TINT, zorder=0, linewidth=0)
    for gx in (-4, 4, 8):
        ax.axvline(gx, color=GRID, lw=0.6, zorder=1)
    ax.axvline(0.0, color=SLATE, lw=0.9, zorder=2)

    for j, seed in enumerate(seeds):
        off = 0.175 if j == 0 else -0.175
        s = SEED_STYLE[seed]
        vals = [deg["deviation_percent"].get(key, {}).get(seed) for key, _ in NODE_ROWS]
        xs = [0.0 if v is None else v for v in vals]
        ax.hlines(
            y=y + off,
            xmin=0.0,
            xmax=xs,
            color=s["color"],
            lw=1.5,
            zorder=3,
            capstyle="round",
        )
        ax.scatter(
            xs,
            y + off,
            s=MARKER_SIZE,
            color=s["color"],
            marker=s["marker"],
            edgecolors="white",
            linewidths=0.5,
            zorder=5,
            label=s["label"],
        )
        for yi, v in zip(y, vals):
            if v is None:
                continue
            ax.text(
                v + (0.45 if v >= 0 else -0.45),
                yi + off,
                _signed_pct(v),
                ha="left" if v >= 0 else "right",
                va="center",
                fontsize=8,
                color=s["color"],
                zorder=6,
            )

    ax.set_yticks(y, [label for _, label in NODE_ROWS], fontsize=8.5)
    ax.set_ylim(-0.62, len(NODE_ROWS) - 0.5)
    ax.set_xlim(x_lo, x_hi)
    ax.set_xticks([-4, 0, 4, 8])
    ax.set_xticklabels(["\u22124", "0", "4", "8"])
    ax.set_xlabel("Deviation (uniform \u2212 gamma) / gamma (%)", labelpad=4)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.spines["left"].set_visible(False)
    ax.invert_yaxis()
    _header_legend(ax, ncol=2, columnspacing=2.0, handlelength=0.95, handletextpad=0.35)


def main():
    tier1 = load(TIER1_JSON)
    deg = load(DEG_JSON)

    # 手工布局：两栏轴位置完全一致，页眉（图例行 + 标题行）因此天然对齐
    fig = plt.figure(figsize=(6.8, 2.78))
    ax_a = fig.add_axes([0.082, 0.205, 0.395, 0.545])
    ax_b = fig.add_axes([0.655, 0.205, 0.330, 0.545])
    panel_a(ax_a, tier1)
    panel_b(ax_b, deg)

    for ax, title in (
        (ax_a, "(a) Wrapper agreement with native runs"),
        (ax_b, "(b) Gamma-to-uniform degradation"),
    ):
        ax.text(
            0.0,
            1.230,
            title,
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            color=INK,
            ha="left",
            va="bottom",
        )

    # 仅保存图像，不做任何屏幕预览
    fig.savefig(f"{SCRIPT_NAME}.pdf", bbox_inches="tight", pad_inches=0.02)
    png = f"{SCRIPT_NAME}.png"
    fig.savefig(png, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    from PIL import Image  # 合成到不透明白底，写入 300 dpi 元数据

    with Image.open(png) as im:
        if im.mode != "RGB":
            bg = Image.new("RGB", im.size, (255, 255, 255))
            bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
            bg.save(png, dpi=(600, 600))
    plt.close(fig)
    print(f"saved: {SCRIPT_NAME}.pdf / {SCRIPT_NAME}.png")


if __name__ == "__main__":
    main()
