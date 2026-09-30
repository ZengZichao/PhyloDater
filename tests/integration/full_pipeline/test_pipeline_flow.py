#!/usr/bin/env python3
"""
PhyloDater 全流程定年测试
使用PhyloDater Pipeline，每个工具独立文件夹，完整文件输出
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

# 项目根目录
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from phylodater.core.pipeline import DatingPipeline, PipelineConfig
from phylodater.infrastructure.configuration import ToolConfig
from phylodater.infrastructure.safe_io import safe_writer
from phylodater.models import CalibrationPoint, FixedAgeConstraint, PhylogeneticTree

# 输出目录（项目内 output/integration_pipeline，避免写入桌面）
BASE_OUTPUT_DIR = (
    PROJECT_ROOT / "output" / "integration_pipeline" / "phylodater_full_test"
)
BASE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 日志配置
LOG_FILE = BASE_OUTPUT_DIR / "phylodater_main.log"

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(LOG_FILE, mode="w"), logging.StreamHandler()],
)
logger = logging.getLogger("phylodater_test")

# 测试数据
NEWICK_TREE = "((human:0.03,(chimp:0.01,bonobo:0.01):0.02):0.05,(gorilla:0.06,(orangutan:0.04,sumatran:0.04):0.02):0.02):0.02;"

ALIGNMENT_CONTENT = """>human
AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTAT
>chimp
TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA
>bonobo
AAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTACAATGGCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAA
>gorilla
ACAGCTCAGGAGCCAGTCCCCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGGTATCGATACGTAGGAGGTGTGCCTTCGTAGGCTGTTTCTC
>orangutan
AGGACGCCCAACTATTCTTTCCAATCCTACATCTGTTTCTTGCGTCGTAGCGGGACCCTCCATTGTTACTTATTAGGTTCTCGTTATGTCTCATAATCTC
>sumatran
AGTGCTGGTGTGATAAGCAAACCACCCTACTGGCACGAAGTTCACAGAAGTGAGATTATGTCTCGTTTGGCAGTCTTGATGCTCGGGGGACACTTCTTTA
"""


def create_test_files():
    """创建测试文件"""
    data_dir = BASE_OUTPUT_DIR / "input_data"
    data_dir.mkdir(exist_ok=True)

    tree_file = data_dir / "tree.nwk"
    with safe_writer(tree_file) as f:
        f.write(NEWICK_TREE)
    logger.info(f"Created tree file: {tree_file}")

    alignment_file = data_dir / "alignment.fasta"
    with safe_writer(alignment_file) as f:
        f.write(ALIGNMENT_CONTENT)
    logger.info(f"Created alignment file: {alignment_file}")

    # 创建校准点 - 设置resolved_taxa以跳过自动解析
    # 注意：根节点校准不放入calibrations列表，通过ToolConfig.mcmctree.root_age设置
    # MCMCTree 4.10.8+不支持在树文件中添加根节点校准，必须通过控制文件RootAge参数
    calibrations = [
        CalibrationPoint(
            name="HumanChimp",
            mrca_leaf_pair=("human", "chimp"),
            resolved_taxa=["human", "chimp"],  # 预解析
            age_constraint=FixedAgeConstraint(fixed_age=6.0),
        ),
    ]

    logger.info(f"Created {len(calibrations)} calibration points")
    for cal in calibrations:
        logger.info(f"  - {cal.name}: {cal.age_constraint}")

    return tree_file, alignment_file, calibrations


def plot_tree_from_result(result, method_name, output_path, root_age=15.0):
    """从DatingResult绘制时间树"""
    import re

    from ete3 import Tree

    newick = result.dated_tree_newick
    clean_newick = re.sub(r"\[.*?\]", "", newick)

    ete_tree = Tree(clean_newick, format=1)

    max_dist = 0
    for leaf in ete_tree.get_leaves():
        dist = ete_tree.get_distance(leaf)
        if dist > max_dist:
            max_dist = dist

    if max_dist == 0:
        max_dist = 1.0

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
            child_ys = []
            for child in node.children:
                child_y = assign_coords(child)
                child_ys.append(child_y)
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
    ax.set_ylim(0, 7)
    ax.invert_xaxis()
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_xlabel("Age (Ma)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Taxa", fontsize=12, fontweight="bold")
    ax.set_title(
        f"{method_name} Dating Result\n(Root Age: {root_age:.1f} Ma)",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(loc="upper left", fontsize=9)

    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    logger.info(f"Saved tree plot: {output_path}")


def run_single_method(method_name, tree_file, alignment_file, calibrations):
    """运行单个定年方法"""
    logger.info("=" * 60)
    logger.info(f"Running method: {method_name}")
    logger.info("=" * 60)

    method_dir = BASE_OUTPUT_DIR / method_name
    method_dir.mkdir(exist_ok=True)

    method_log_file = method_dir / f"{method_name}.log"
    method_logger = logging.getLogger(f"phylodater.{method_name}")
    method_logger.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(method_log_file, mode="w")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    method_logger.addHandler(file_handler)

    try:
        with open(tree_file) as f:
            tree_newick = f.read().strip()
        tree = PhylogeneticTree(_canonical_newick=tree_newick)

        # 创建ToolConfig并设置MCMCTree的root_age
        tool_config = ToolConfig()
        tool_config.mcmctree.root_age = (
            20.0  # 根节点年龄上界（与MaximumAgeConstraint一致）
        )

        config = PipelineConfig(
            methods=[method_name],
            output_dir=method_dir,
            threads=1,
            skip_failed=True,
            seed=42,
        )

        pipeline = DatingPipeline(config, tool_config=tool_config)

        method_logger.info(f"Starting {method_name} analysis")
        method_logger.info(f"Tree: {tree_newick}")
        method_logger.info(f"Alignment: {alignment_file}")
        method_logger.info(f"Calibrations: {len(calibrations)}")

        results = pipeline.run(tree, alignment_file, calibrations)

        if method_name in results:
            result = results[method_name]

            result_tree_file = method_dir / f"{method_name}_result.nwk"
            with safe_writer(result_tree_file) as f:
                f.write(result.dated_tree_newick)
            method_logger.info(f"Saved result tree: {result_tree_file}")

            ages_file = method_dir / f"{method_name}_ages.txt"
            with safe_writer(ages_file) as f:
                f.write(f"Node Ages ({method_name})\n")
                f.write("=" * 40 + "\n")
                for node_name, age_estimate in result.node_ages.items():
                    f.write(f"{node_name}: {age_estimate.mean_age:.2f} Ma")
                    if age_estimate.ci_lower and age_estimate.ci_upper:
                        f.write(
                            f" (95% CI: {age_estimate.ci_lower:.2f}-{age_estimate.ci_upper:.2f})"
                        )
                    f.write("\n")
            method_logger.info(f"Saved node ages: {ages_file}")

            plot_file = method_dir / f"{method_name}_tree.png"
            plot_tree_from_result(result, method_name, plot_file)

            method_logger.info(f"{method_name} completed successfully")
            logger.info(f"✓ {method_name} completed")
            return True
        else:
            method_logger.error(f"{method_name} failed - no result returned")
            logger.error(f"✗ {method_name} failed")
            return False

    except Exception as e:
        method_logger.error(f"{method_name} failed with exception: {e}", exc_info=True)
        logger.error(f"✗ {method_name} failed: {e}")
        return False
    finally:
        method_logger.removeHandler(file_handler)


def main():
    """运行所有测试"""
    logger.info("=" * 60)
    logger.info("PhyloDater Full Pipeline Test")
    logger.info("=" * 60)
    logger.info(f"Output directory: {BASE_OUTPUT_DIR}")
    logger.info(f"Main log file: {LOG_FILE}")
    logger.info(f'Test time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

    tree_file, alignment_file, calibrations = create_test_files()

    methods = ["pathd8", "mcmctree", "wlogdate", "mdcat", "pyr8s", "lsd2"]

    results = {}

    for method in methods:
        success = run_single_method(method, tree_file, alignment_file, calibrations)
        results[method] = success

    logger.info("")
    logger.info("=" * 60)
    logger.info("Test Results Summary")
    logger.info("=" * 60)

    for method, success in results.items():
        status = "✓ Success" if success else "✗ Failed"
        logger.info(f"  {method:15s}: {status}")

    success_count = sum(1 for v in results.values() if v)
    total_count = len(results)

    logger.info("")
    logger.info(f"Total: {success_count}/{total_count} successful")
    logger.info(f"Results saved in: {BASE_OUTPUT_DIR}")
    logger.info(f"Main log: {LOG_FILE}")

    logger.info("")
    logger.info("Generated files:")
    for method_dir in sorted(BASE_OUTPUT_DIR.iterdir()):
        if method_dir.is_dir() and method_dir.name != "input_data":
            logger.info(f"  {method_dir.name}/")
            for file in sorted(method_dir.iterdir()):
                logger.info(f"    - {file.name} ({file.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
