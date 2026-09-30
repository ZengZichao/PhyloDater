"""
PhyloDater 全流程定年测试

调用所有外部定年工具并绘制结果图片。
需要安装 PATHd8、MCMCTree、wLogDate、MD-Cat、pyr8s、IQ-TREE2。
"""

import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import LineCollection


def _wlogdate_available() -> bool:
    """检查 wLogDate 的 Python 依赖 logdate 是否可导入。"""
    try:
        import logdate  # noqa: F401

        return True
    except ImportError:
        return False


# 测试树
NEWICK_TREE = "((human:0.03,(chimp:0.01,bonobo:0.01):0.02):0.05,(gorilla:0.06,(orangutan:0.04,sumatran:0.04):0.02):0.02):0.02;"
NEWICK_TREE_SHORT = "((hum:0.03,(chi:0.01,bon:0.01):0.02):0.05,(gor:0.06,(ora:0.04,sum:0.04):0.02):0.02):0.02;"

TAXA = ["human", "chimp", "bonobo", "gorilla", "orangutan", "sumatran"]
TAXA_SHORT = ["hum", "chi", "bon", "gor", "ora", "sum"]


def generate_sequences(length=500, short_names=False):
    """生成随机 DNA 序列"""
    import random

    random.seed(42)
    bases = ["A", "C", "G", "T"]
    taxa = TAXA_SHORT if short_names else TAXA
    return {
        taxon: "".join(random.choice(bases) for _ in range(length)) for taxon in taxa
    }


def plot_tree(newick, title, output_path, root_age=15.0, method_name=""):
    """绘制时间树"""
    import re

    from ete3 import Tree

    clean_newick = re.sub(r"\[.*?\]", "", newick)
    ete_tree = Tree(clean_newick, format=1)

    max_dist = max(ete_tree.get_distance(leaf) for leaf in ete_tree.get_leaves()) or 1.0

    coords = {}
    y_counter = [0]

    def get_age(node):
        if node.is_leaf():
            return 0.0
        dist = ete_tree.get_distance(node)
        return root_age * (1 - dist / max_dist)

    def assign_coords(node):
        if node.is_leaf():
            y_counter[0] += 1
            y = float(y_counter[0])
        else:
            child_ys = [assign_coords(child) for child in node.children]
            y = np.mean(child_ys)
        x = get_age(node)
        coords[node] = (x, y)
        return y

    assign_coords(ete_tree)

    segments = []
    for node in ete_tree.traverse():
        if node.is_root():
            continue
        if node not in coords or node.up not in coords:
            continue
        x_child, y_child = coords[node]
        x_parent, _ = coords[node.up]
        segments.append([(x_parent, y_child), (x_child, y_child)])

    for node in ete_tree.traverse():
        if node.is_leaf():
            continue
        if node not in coords:
            continue
        child_ys = [coords[child][1] for child in node.children if child in coords]
        if len(child_ys) >= 2:
            x = coords[node][0]
            segments.append([(x, min(child_ys)), (x, max(child_ys))])

    fig, ax = plt.subplots(figsize=(14, 9))
    lc = LineCollection(segments, colors="#2C3E50", linewidths=2.5, zorder=1)
    ax.add_collection(lc)

    leaf_plotted = False
    internal_plotted = False
    for node, (x, y) in coords.items():
        if node.is_leaf():
            label = "Leaf" if not leaf_plotted else None
            ax.scatter(
                x,
                y,
                s=100,
                c="#E74C3C",
                zorder=3,
                edgecolors="white",
                linewidths=1.5,
                label=label,
            )
            ax.annotate(
                node.name,
                xy=(x, y),
                xytext=(8, 0),
                textcoords="offset points",
                fontsize=10,
                fontweight="bold",
                fontstyle="italic",
            )
            leaf_plotted = True
        else:
            label = "Internal" if not internal_plotted else None
            ax.scatter(
                x,
                y,
                s=60,
                c="#3498DB",
                zorder=3,
                edgecolors="white",
                linewidths=1,
                label=label,
            )
            internal_plotted = True

    ax.set_xlim(-root_age * 0.05, root_age * 1.1)
    ax.set_ylim(0, len(TAXA) + 1)
    ax.invert_xaxis()
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_xlabel("Age (Ma)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Taxa", fontsize=12, fontweight="bold")
    ax.set_title(
        f"{title}\n(Method: {method_name}, Root Age: {root_age:.1f} Ma)",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(loc="upper left", fontsize=9)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()


def test_pathd8(tmp_path):
    """测试 PATHd8"""
    if not shutil.which("PATHd8"):
        pytest.skip("PATHd8 not installed")

    tmpdir = tmp_path / "pathd8"
    tmpdir.mkdir()
    infile = tmpdir / "test.infile"
    outfile = tmpdir / "test.outfile"

    infile.write_text(
        "Sequence length = 500;\n"
        + NEWICK_TREE
        + "\n\n"
        + "mrca: human, chimp, fixage=6;\n"
        + "mrca: human, gorilla, minage=8;\n"
        + "mrca: human, gorilla, maxage=10;\n"
        + "name of mrca: human, chimp, name=HumanChimp;\n"
        + "name of mrca: human, gorilla, name=GreatApes;\n"
    )

    result = subprocess.run(
        ["PATHd8", str(infile), str(outfile)],
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert outfile.exists(), f"PATHd8 did not create output: {result.stderr}"
    content = outfile.read_text()
    tree_line = None
    for line in content.split("\n"):
        if "d8 tree" in line and ":" in line:
            tree_line = line.split(":", 1)[1].strip()
            break

    assert tree_line and tree_line.endswith(";"), "PATHd8 output tree not found"
    plot_tree(
        tree_line,
        "PATHd8 Dating Result",
        tmpdir / "pathd8.png",
        root_age=15,
        method_name="PATHd8",
    )


def test_mcmctree(tmp_path):
    """测试 MCMCTree。未安装时自动跳过，安装后在本地执行轻量参数版本。"""
    if not shutil.which("mcmctree"):
        pytest.skip("mcmctree not installed")

    tmpdir = tmp_path / "mcmctree"
    tmpdir.mkdir()
    seqs = generate_sequences(100, short_names=True)

    seq_file = tmpdir / "mcmctree.phy"
    seq_file.write_text(
        f"  {len(TAXA_SHORT)} 100\n"
        + "".join(f"{name}   {seq}\n" for name, seq in seqs.items())
    )

    tree_file = tmpdir / "mcmctree.tree"
    tree_file.write_text(f"  {len(TAXA_SHORT)}  1\n{NEWICK_TREE_SHORT}")

    ctl_file = tmpdir / "mcmctree.ctl"
    # PAML >= 4.10 requires the birth-death construction flag ('c' conditional /
    # 'm' multiplicative) as the 4th BDparas field; without it mcmctree aborts
    # with "BDparas: expect flag for birth-death process prior". This mirrors the
    # control file produced by phylodater.adapters.mcmctree_method.
    ctl_file.write_text(
        "seed = -1\nseqfile = mcmctree.phy\ntreefile = mcmctree.tree\n"
        "outfile = mcmctree.out\nndata = 1\nseqtype = 0\nusedata = 1\nclock = 2\n"
        "RootAge = '<15.0'\nmodel = 0\nalpha = 0\nncatG = 5\ncleandata = 0\n"
        "BDparas = 1 1 0.1 c\nkappa_gamma = 6 2\nalpha_gamma = 1 1\n"
        "rgene_gamma = 2 20\nsigma2_gamma = 1 10\n"
        "print = 0\nburnin = 100\nsampfreq = 2\nnsample = 200\n"
    )

    result = subprocess.run(
        ["mcmctree", "mcmctree.ctl"],
        capture_output=True,
        text=True,
        timeout=300,
        cwd=str(tmpdir),
    )

    out_file = tmpdir / "mcmctree.out"
    assert out_file.exists(), f"MCMCTree did not create output: {result.stderr}"
    content = out_file.read_text()
    assert "Species tree" in content or len(content) > 1000

    plot_tree(
        NEWICK_TREE,
        "MCMCTree Dating Result",
        tmpdir / "mcmctree.png",
        root_age=15,
        method_name="MCMCTree",
    )


def test_wlogdate(tmp_path):
    """测试 wLogDate"""
    wlogdate_script = shutil.which("launch_wLogDate.py")
    if not wlogdate_script or not _wlogdate_available():
        pytest.skip(
            "wLogDate not installed or its Python dependency 'logdate' is missing"
        )

    tmpdir = tmp_path / "wlogdate"
    tmpdir.mkdir()
    tree_file = tmpdir / "input.tre"
    tree_file.write_text(NEWICK_TREE)
    output_file = tmpdir / "output.tre"

    result = subprocess.run(
        ["python", wlogdate_script, "-i", str(tree_file), "-b", "-o", str(output_file)],
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert output_file.exists(), f"wLogDate did not create output: {result.stderr}"
    content = output_file.read_text()
    assert "t=" in content
    plot_tree(
        content.strip(),
        "wLogDate Dating Result",
        tmpdir / "wlogdate.png",
        root_age=15,
        method_name="wLogDate",
    )


def test_mdcat(tmp_path):
    """测试 MD-Cat"""
    mdcat_script = shutil.which("md_cat.py")
    if not mdcat_script:
        pytest.skip("MD-Cat not installed")

    tmpdir = tmp_path / "mdcat"
    tmpdir.mkdir()
    tree_file = tmpdir / "input.tre"
    tree_file.write_text(NEWICK_TREE)
    output_file = tmpdir / "output.tre"

    result = subprocess.run(
        [
            "python",
            mdcat_script,
            "-i",
            str(tree_file),
            "-b",
            "-o",
            str(output_file),
            "-p",
            "10",
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert output_file.exists(), f"MD-Cat did not create output: {result.stderr}"
    content = output_file.read_text()
    assert content.strip()
    plot_tree(
        content.strip(),
        "MD-Cat Dating Result",
        tmpdir / "mdcat.png",
        root_age=15,
        method_name="MD-Cat",
    )


def test_pyr8s(tmp_path):
    """测试 pyr8s"""
    if not shutil.which("pyr8s"):
        pytest.skip("pyr8s not installed")

    tmpdir = tmp_path / "pyr8s"
    tmpdir.mkdir()
    nexus_file = tmpdir / "test.nex"
    nexus_file.write_text(
        "#NEXUS\n\n"
        "begin trees;\n"
        f"tree test = [&U] {NEWICK_TREE}\n"
        "end;\n\n"
        "begin rates;\n\n"
        "blformat nsites=500 lengths=persite;\n\n"
        "collapse;\n\n"
        "mrca HumanChimp human chimp;\n"
        "mrca GreatApes human gorilla;\n\n"
        "fixage taxon=HumanChimp age=6;\n\n"
        "divtime method=np algorithm=pl;\n"
        "showage;\n"
        "describe plot=chronogram;\n\n"
        "end;\n"
    )

    result = subprocess.run(
        ["pyr8s", str(nexus_file)], capture_output=True, text=True, timeout=60
    )

    assert (
        "CHRONOGRAM" in result.stdout or "SHOWAGE" in result.stdout
    ), f"pyr8s output missing expected markers: {result.stderr}"
    plot_tree(
        NEWICK_TREE,
        "pyr8s Dating Result",
        tmpdir / "pyr8s.png",
        root_age=15,
        method_name="pyr8s",
    )


def test_lsd2(tmp_path):
    """测试 LSD2 (IQ-TREE)"""
    if not shutil.which("iqtree"):
        pytest.skip("IQ-TREE2 not installed")

    tmpdir = tmp_path / "lsd2"
    tmpdir.mkdir()
    seqs = generate_sequences(500)

    seq_file = tmpdir / "test.phy"
    seq_file.write_text(
        f"{len(TAXA)} 500\n"
        + "".join(f"{name:10s}{seq}\n" for name, seq in seqs.items())
    )

    tree_file = tmpdir / "test.tree"
    tree_file.write_text(NEWICK_TREE)

    date_file = tmpdir / "test.date"
    date_file.write_text(
        "human 2020\nchimp 2019\nbonobo 2018\n"
        "gorilla 2015\norangutan 2014\nsumatran 2013\n"
    )

    result = subprocess.run(
        [
            "iqtree",
            "-s",
            str(seq_file),
            "-t",
            str(tree_file),
            "--date",
            str(date_file),
            "-m",
            "JC69",
            "--prefix",
            str(tmpdir / "lsd2"),
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )

    nwk_file = tmpdir / "lsd2.timetree.nwk"
    assert nwk_file.exists(), f"LSD2 did not create timetree: {result.stderr}"
    content = nwk_file.read_text()
    assert content.strip()
    plot_tree(
        content.strip(),
        "LSD2 (IQ-TREE) Dating Result",
        tmpdir / "lsd2.png",
        root_age=15,
        method_name="LSD2",
    )
