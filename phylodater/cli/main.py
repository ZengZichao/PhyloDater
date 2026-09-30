"""
PhyloDater 命令行接口

子命令架构：
  phylodater dating   - 执行定年分析
  phylodater check    - 系统自检
  phylodater fold     - 折叠树（按分类群聚合）〔实验性：尚未实现，调用即报错〕

主解析器仅包含通用参数，各子命令有独立的参数集。
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional

from .. import show_banner
from ..core import DatingMethodRegistry
from ..core.exceptions import (
    AlignmentValidationError,
    CalibrationError,
    ConfigurationError,
    ConsistencyError,
    NegativeBranchLengthError,
    TreeValidationError,
    UnknownMethodError,
)
from ..core.pipeline import DatingPipeline, ParallelPipeline, PipelineConfig
from ..infrastructure import Configuration, LogLevel, get_logger

if TYPE_CHECKING:  # 只供标注使用，避免 CLI 启动期多一次导入
    from ..infrastructure.logging import PhyloDaterLogger
from ..models import PhylogeneticTree
from ..services import CalibrationLoader
from ..services.auto_calibrator import AutoCalibrator, parse_auto_calibrate_args
from ..services.deep_validator import DeepValidator
from ..services.input_validator import InputValidator

# ==================== 退出码定义 ====================
EXIT_SUCCESS = 0  # 成功
EXIT_RUNTIME_ERROR = 1  # 运行错误（依赖缺失、自检失败等）
EXIT_USAGE_ERROR = 2  # 参数或命令行使用错误
EXIT_DATA_ERROR = 3  # 输入数据格式或内容错误
EXIT_INTERRUPTED = 130  # 用户中断

# ==================== 未实现能力的对外表述 ====================
# 审阅项 C-36：`fold` 有完整的参数面，但处理函数只 raise NotImplementedError。
# 这本身是"响亮失败"，问题在于它此前**被当作已交付能力呈现**：出现在 --help
# 的正常列表里、参数被完整接受并回显，用户要跑到最后一行才知道功能不存在。
# 修复放在 CLI 层（不伪造实现、也不删掉那条正确的未实现错误）：
# 把"实验性/未实现"这一事实前移到用户在**决定调用之前**就能看到的位置
# —— 主 help 的子命令清单、`fold --help` 的 help/description，以及执行时
# 先于参数回显的第一条日志。措辞集中在此常量里，避免各处漂移。
EXPERIMENTAL_NOT_IMPLEMENTED = "实验性，尚未实现 (experimental, not yet implemented)"


def _get_version_string() -> str:
    """获取详细的版本信息，包括提交哈希和依赖版本"""
    from .. import __license__, __update_date__, __version__

    version_parts = [
        f"PhyloDater {__version__}",
        f"License: {__license__}",
        f"Updated: {__update_date__}",
    ]

    # 静态版本不含 ``+``；仅当版本串带 local 段（如自定义构建）时才解析其中的 git 哈希。
    if "+" in __version__:
        git_part = __version__.split("+")[-1]
        if git_part.startswith("g") and len(git_part) >= 7:
            version_parts.append(f"git: {git_part[1:]}")
    else:
        # 回退：在 Git 工作区内用 git 命令读取提交短哈希（发布版经 pip 安装时通常无 git，静默跳过）。
        try:
            result = subprocess.run(
                ["git", "rev-parse", "--short", "HEAD"],
                capture_output=True,
                text=True,
                timeout=5,
                cwd=Path(__file__).parent.parent.parent,
            )
            if result.returncode == 0:
                version_parts.append(f"git: {result.stdout.strip()}")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

    # 获取关键依赖版本
    deps = []
    for mod_name, label in [
        ("Bio", "BioPython"),
        ("ete3", "ETE3"),
        ("dendropy", "DendroPy"),
        ("numpy", "NumPy"),
    ]:
        try:
            mod = __import__(mod_name)
            deps.append(f"{label}={getattr(mod, '__version__', '?')}")
        except (ImportError, AttributeError):
            pass

    if deps:
        version_parts.append(f"deps: {', '.join(deps)}")

    return "\n".join(version_parts)


def _validate_existing_path(value: str) -> Path:
    """验证路径是否存在，不存在时以 EXIT_DATA_ERROR 退出"""
    path = Path(value)
    if not path.exists():
        sys.stderr.write(f"错误: 路径不存在: {path}\n")
        sys.exit(EXIT_DATA_ERROR)
    return path


def _validate_config_path(path_str: str) -> Path:
    """验证配置文件路径，不存在时以 EXIT_USAGE_ERROR 退出"""
    p = Path(path_str)
    if not p.exists():
        print(f"错误: 配置文件不存在: {path_str}", file=sys.stderr)
        sys.exit(EXIT_USAGE_ERROR)
    return p


def _validate_positive_int(value: str) -> int:
    """验证正整数参数"""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"必须是整数，得到: {value}")
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"必须是正整数，得到: {value}")
    return ivalue


def _clamp_threads(value: int) -> int:
    """将线程数钳制到 [1, os.cpu_count()] 区间。

    ``--threads``/``-j`` 的负值/非整数已由 :func:`_validate_positive_int` 在解析阶段
    拒绝；此处仅处理超大值，避免传入远超物理核数导致 CPU 过度订阅。
    """
    max_threads = os.cpu_count() or 1
    return max(1, min(int(value), max_threads))


def _parse_method_args(raw: str) -> Dict[str, str]:
    """
    解析 --method-args 键值对

    格式: "key1=value1,key2=value2" 或 "key1=value1 key2=value2"
    """
    result: Dict[str, str] = {}
    if not raw:
        return result

    # 支持逗号或空格分隔
    for pair in raw.replace(" ", ",").split(","):
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            raise argparse.ArgumentTypeError(
                f"method-args 格式错误: '{pair}'，应为 key=value 格式"
            )
        k, v = pair.split("=", 1)
        result[k.strip()] = v.strip()

    return result


# ==================== 子命令创建 ====================


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """添加通用参数（所有子命令共享）"""
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="日志级别 (默认: INFO)",
    )

    parser.add_argument(
        "--log-file", type=Path, default=None, help="日志文件路径（UTF-8，实时 flush）"
    )

    parser.add_argument(
        "--taxonomy-file",
        type=Path,
        default=None,
        help="外部分类映射文件 (TSV/CSV: name, taxonomy_string)",
    )

    parser.add_argument(
        "--taxonomy-levels",
        type=str,
        nargs="+",
        default=None,
        metavar="PREFIX:RANK",
        help='自定义分类级别前缀 (如 "k:kingdom ss:subspecies")',
    )

    parser.add_argument(
        "--taxonomy-delimiter-mode",
        choices=["reverse", "greedy", "segment"],
        default="reverse",
        help="格式A解析策略 (默认: reverse)",
    )

    parser.add_argument(
        "--taxonomy-source-priority",
        choices=["embedded", "table"],
        default="table",
        help="分类学来源优先级 (默认: table)",
    )

    parser.add_argument(
        "--table-sep", type=str, default=None, help="强制指定分类学表格分隔符"
    )

    parser.add_argument(
        "--taxonomy-table-sep",
        default=";",
        help="分类学字符串内部的级别分隔符（默认 ;，可设为 | 等）",
    )

    parser.add_argument(
        "--ignore-malformed",
        action="store_true",
        default=False,
        help="遇到畸形输入时跳过（默认终止）",
    )

    parser.add_argument("--version", action="version", version=_get_version_string())


def _create_dating_subparser(
    subparsers: argparse._SubParsersAction, common_parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """创建 dating 子命令"""
    epilog = """示例用法:
  # 使用 MCMCTree 定年
  phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o output/ --method mcmctree

  # 使用 treePL 定年，指定平滑参数
  phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o output/ --method treepl --method-args "smooth=100"

  # 使用 MCMCTree + PAML 4.8
  phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml \\
      -o output/ --method mcmctree --method-args "paml_version=4.8"

  # 使用 r8s 定年
  phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o output/ --method r8s

  # 干运行（仅验证）
  phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o output/ --method mcmctree --dry-run

  # 自动添加校准点
  phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml \\
      -o output/ --method mcmctree --auto-calibrate "Cyanobacteriota:2500"

方法专有参数 (--method-args):
  mcmctree:  clock, num_runs, burnin, nsample, sampfreq,
             paml_version, paml_path, mcmctree_bin,
             psrf_threshold, ess_threshold
  treepl:    initial_smooth, cvstart, cvstop, cviter, cvmultstep
  pathd8:    (无额外参数；任何取值都会被拒绝)
  r8s:       method(LF/NPRS/PL), smoothing, algorithm, r8s_bin
  lsd2:      model, iqtree_bin, date_ci, clock_sd, date_options, date_outlier
  wlogdate:  backward_time, num_replicates, root_time, leaf_time,
             sequence_length, max_iter, pseudocount, zero_branch_len
  mdcat:     ncat, nrep, max_iter, backward_time, as_date,
             use_direct_import, ci_nboots, ci_plower, ci_pupper, annotate_level

未知参数一律以退出码 3 拒绝（不静默丢弃），见 Configuration.set_cli_override。
"""

    dating_parser: argparse.ArgumentParser = subparsers.add_parser(
        "dating",
        help="执行定年分析",
        description="使用指定的定年方法对系统发育树进行分歧时间估算",
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[common_parser],
    )

    # === 输入文件 ===
    input_group = dating_parser.add_argument_group("输入文件")

    input_group.add_argument(
        "--tree",
        "-t",
        required=True,
        type=_validate_existing_path,
        help="输入系统发育树文件 (Newick/Nexus 格式)",
    )

    input_group.add_argument(
        "--sequence",
        "-s",
        required=True,
        type=_validate_existing_path,
        help="输入序列比对文件 (FASTA 格式)",
    )

    input_group.add_argument(
        "--calibrations",
        "-c",
        required=True,
        type=_validate_existing_path,
        help="校准点配置文件 (YAML 格式)",
    )

    input_group.add_argument(
        "--output", "-o", required=True, type=Path, help="输出目录"
    )

    input_group.add_argument(
        "--config", type=_validate_config_path, help="配置文件路径 (YAML 格式)"
    )

    input_group.add_argument(
        "--generate-config",
        type=Path,
        metavar="OUTPUT_PATH",
        help="生成默认配置文件到指定路径并退出",
    )

    # === 定年方法 ===
    method_group = dating_parser.add_argument_group("定年方法")

    method_group.add_argument(
        "--method",
        required=True,
        choices=DatingMethodRegistry.list_methods(),
        help="定年方法名称",
    )

    method_group.add_argument(
        "--method-args",
        type=_parse_method_args,
        default={},
        metavar="KEY=VALUE,...",
        help='方法专有参数，逗号分隔的键值对 (如 "clock=2,num_runs=2")',
    )

    # === 树处理 ===
    tree_group = dating_parser.add_argument_group("树处理")

    tree_group.add_argument(
        "--reroot",
        choices=["outgroup", "midpoint", "mad", "none"],
        default="none",
        help="定根策略 (默认: none)",
    )

    tree_group.add_argument("--outgroup", help="外群名称（用于 outgroup 定根）")

    tree_group.add_argument(
        "--multi-tree-mode",
        choices=["split", "first", "last", "random", "ask", "error"],
        default="error",
        help="多棵树处理模式 (默认: error, ask=交互选择后退出码2)",
    )

    # === 校准 ===
    cal_group = dating_parser.add_argument_group("校准参数")

    cal_group.add_argument(
        "--auto-calibrate",
        nargs="+",
        metavar="TAXON:AGE",
        default=None,
        help='自动添加校准点 (如 "Cyanobacteriota:2500")',
    )

    # === 执行控制 ===
    exec_group = dating_parser.add_argument_group("执行控制")

    exec_group.add_argument(
        "--threads",
        "-j",
        type=_validate_positive_int,
        default=4,
        help="并行线程数 (默认: 4)",
    )

    exec_group.add_argument("--seed", type=int, default=None, help="全局随机种子")

    exec_group.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="干运行模式：仅验证配置和依赖",
    )

    exec_group.add_argument(
        "--mol-type",
        choices=["DNA", "RNA", "protein"],
        default=None,
        help="序列分子类型（默认自动检测）",
    )

    exec_group.add_argument(
        "--no-cross-check",
        action="store_true",
        default=False,
        help="关闭树与序列的末端标签交叉验证",
    )

    exec_group.add_argument(
        "--skip-length-check",
        action="store_true",
        default=False,
        help="跳过负分支长度检查和序列长度一致性检查",
    )

    exec_group.add_argument(
        "--low-memory", action="store_true", default=False, help="低内存模式"
    )

    exec_group.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help=(
            "从检查点恢复执行。已缓存结果只有在**输入指纹一致**时才会被复用："
            "树、比对、校准点、生效配置或软件路径任一变化，或检查点里没有可核验"
            "的指纹，都会全量重跑而非照搬上次结果 (审阅项 B-23)"
        ),
    )

    exec_group.add_argument(
        "--clear-checkpoint",
        action="store_true",
        default=False,
        help="清除检查点文件后重新执行",
    )

    # === 输出控制 ===
    out_group = dating_parser.add_argument_group("输出控制")

    out_group.add_argument(
        "--force", action="store_true", default=False, help="覆盖已存在的输出文件"
    )

    out_group.add_argument(
        "--no-clobber", action="store_true", default=False, help="跳过已存在的输出文件"
    )

    out_group.add_argument(
        "--strip-annotations",
        action="store_true",
        default=False,
        help="输出时丢弃树注释（Bootstrap/NHX）以减小文件体积",
    )

    return dating_parser


def _create_check_subparser(
    subparsers: argparse._SubParsersAction, common_parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """创建 check 子命令"""
    check_parser: argparse.ArgumentParser = subparsers.add_parser(
        "check",
        aliases=["self-test"],
        help="系统自检与依赖验证",
        description="执行系统自检：依赖库检查、示例解析测试、单系群判定验证等",
        parents=[common_parser],
    )
    return check_parser


def _create_fold_subparser(
    subparsers: argparse._SubParsersAction, common_parser: argparse.ArgumentParser
) -> argparse.ArgumentParser:
    """创建 fold 子命令（实验性，尚未实现）。

    保留参数面与 ``help``（便于用户读到明确声明、也便于实现后直接接上），
    但在 ``help``/``description``/``epilog`` 三处都标注"尚未实现"，让"这项
    能力当前不可用"出现在用户决定调用**之前**（审阅项 C-36）。
    """
    fold_parser: argparse.ArgumentParser = subparsers.add_parser(
        "fold",
        help=f"按分类群折叠树（聚合末端节点）【{EXPERIMENTAL_NOT_IMPLEMENTED}】",
        description=(
            "将系统发育树按分类学级别折叠，减少末端数量。\n"
            f"注意：本子命令目前【{EXPERIMENTAL_NOT_IMPLEMENTED}】，"
            "参数可解析但不会产出任何折叠树。"
        ),
        parents=[common_parser],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            f"{EXPERIMENTAL_NOT_IMPLEMENTED}\n"
            "调用该子命令会以运行错误码 (exit 1) 结束，不会写出任何折叠树。\n"
            "本功能将在后续版本提供；勿在当前版本中依赖它。"
        ),
    )

    fold_parser.add_argument(
        "--tree", "-t", required=True, type=Path, help="输入系统发育树文件"
    )

    fold_parser.add_argument(
        "--output", "-o", required=True, type=Path, help="输出目录"
    )

    fold_parser.add_argument(
        "--fold-rank",
        choices=["domain", "phylum", "class", "order", "family", "genus"],
        default="family",
        help="折叠到的分类级别 (默认: family)【该子命令尚未实现，此选项目前不生效】",
    )

    return fold_parser


def create_parser() -> argparse.ArgumentParser:
    """创建命令行参数解析器（子命令架构）"""

    # 通用参数单独抽出为 parent parser，使 --log-file/--log-level 等选项
    # 既能放在子命令前，也能放在子命令后。
    common_parser = argparse.ArgumentParser(add_help=False)
    _add_common_arguments(common_parser)

    parser = argparse.ArgumentParser(
        prog="phylodater",
        description="PhyloDater: 多软件并行系统发育定年平台",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[common_parser],
        epilog=f"""子命令:
  dating    执行定年分析
  check     系统自检与依赖验证
  fold      按分类群折叠树【{EXPERIMENTAL_NOT_IMPLEMENTED}，调用即以错误码退出】

使用 'phylodater <subcommand> --help' 查看子命令详细帮助。

依赖库许可证:
  BioPython    - Biopython License (BSD-like)
  ETE3         - GPLv3
  DendroPy     - BSD 3-Clause
  NumPy        - BSD 3-Clause
  Pandas       - BSD 3-Clause
  PyYAML       - MIT License
  Matplotlib   - PSF-based (BSD-compatible)
  psutil       - BSD 3-Clause (可选)
""",
    )

    # 创建子命令解析器
    subparsers = parser.add_subparsers(
        dest="subcommand", title="子命令", help="可用子命令"
    )

    # 注册子命令（继承通用参数）
    _create_dating_subparser(subparsers, common_parser)
    _create_check_subparser(subparsers, common_parser)
    _create_fold_subparser(subparsers, common_parser)

    return parser


# ==================== 子命令执行 ====================


def _setup_logger(
    parsed_args: argparse.Namespace, output_dir: Optional[Path] = None
) -> "PhyloDaterLogger":
    """统一设置日志"""
    log_level = LogLevel(
        parsed_args.log_level if hasattr(parsed_args, "log_level") else "INFO"
    )

    # 两个分支的兼底不同：有 output_dir 时总能落到一个路径，没有时可能真的是 None，
    # 所以这里必须是 Optional[Path]。
    log_file: Optional[Path]
    if output_dir:
        log_file = (
            Path(parsed_args.log_file)
            if hasattr(parsed_args, "log_file") and parsed_args.log_file
            else output_dir / "phylodater.log"
        )
    else:
        log_file = (
            Path(parsed_args.log_file)
            if hasattr(parsed_args, "log_file") and parsed_args.log_file
            else None
        )

    logger = get_logger(log_file)
    logger.level = log_level
    if output_dir:
        logger.set_base_dir(output_dir)
    logger.start()
    return logger


def _run_self_test(logger: "PhyloDaterLogger") -> int:
    """执行自检模式"""
    from ..infrastructure.self_test import SelfTester

    tester = SelfTester()
    all_passed = tester.run_all_checks()
    return EXIT_SUCCESS if all_passed else EXIT_RUNTIME_ERROR


def _run_dating(parsed_args: argparse.Namespace, logger: "PhyloDaterLogger") -> int:
    """执行定年分析子命令"""
    import time

    start_time = time.time()

    # §8.1 - --generate-config: 生成默认配置并退出
    if hasattr(parsed_args, "generate_config") and parsed_args.generate_config:
        config = Configuration()
        config.save_to_file(parsed_args.generate_config)
        logger.info(f"Default configuration saved to {parsed_args.generate_config}")
        return EXIT_SUCCESS

    output_dir = Path(parsed_args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 先记录输出目录再设置日志基准目录，避免路径被相对化显示为 "."
    logger.info(f"输出目录: {output_dir.resolve()}")
    logger.set_base_dir(output_dir)

    # 互斥参数验证
    if parsed_args.force and parsed_args.no_clobber:
        logger.error("--force and --no-clobber are mutually exclusive")
        return EXIT_USAGE_ERROR

    # §14 - 默认检查输出已存在（无--force且无--no-clobber时报ERROR）
    output_marker = output_dir / ".phylodater_run"
    if not parsed_args.force and not parsed_args.no_clobber and output_marker.exists():
        logger.error(
            f"输出目录 {output_dir} 已包含PhyloDater运行结果。"
            "请使用 --force 覆盖或 --no-clobber 跳过已有文件。"
        )
        return EXIT_USAGE_ERROR

    # §11.1 - --clear-checkpoint: 清除检查点
    if hasattr(parsed_args, "clear_checkpoint") and parsed_args.clear_checkpoint:
        from ..infrastructure.checkpoint import PipelineCheckpointManager

        ckpt_mgr = PipelineCheckpointManager(output_dir)
        ckpt_mgr.reset()
        logger.info("Checkpoint cleared")

    logger.section("PhyloDater 定年分析")
    logger.info(f"方法: {parsed_args.method}")
    if parsed_args.method_args:
        logger.info(f"方法参数: {parsed_args.method_args}")

    # 注册信号处理器
    from ..infrastructure.signal_handler import get_shutdown_handler

    shutdown_handler = get_shutdown_handler()
    shutdown_handler.register()
    # 关键修复：绝不能把用户的 output_dir 本体注册为待清理目录——否则 Ctrl+C /
    # atexit 时 shutil.rmtree 会递归删除整个结果目录，造成数据丢失（与
    # ShutdownContext 的修复意图一致：只清理真正独立的临时子目录）。
    # 各方法在 work_dir 下创建的临时目录已通过 pipeline._register_temp_directory
    # 单独注册并清理；此处仅当 output_dir 下确实存在独立临时子目录时才清理它。
    _cli_temp_dir = output_dir / ".phylodater_tmp"
    if _cli_temp_dir.exists():
        shutdown_handler.register_temp_dir(_cli_temp_dir)

    try:
        # 输入验证
        logger.step(1, "输入验证")

        if not _validate_input_files(parsed_args, logger):
            return EXIT_DATA_ERROR

        # 加载树
        logger.step(2, "加载系统发育树")

        if shutdown_handler.is_shutting_down:
            return EXIT_INTERRUPTED

        tree = PhylogeneticTree.from_file(
            parsed_args.tree, multi_tree_mode=parsed_args.multi_tree_mode
        )
        logger.info(f"树包含 {tree.num_tips} 个末端节点")

        # §9.1 - 负分支长度检查 → 退出码3
        skip_length_check = (
            hasattr(parsed_args, "skip_length_check") and parsed_args.skip_length_check
        )
        if not skip_length_check:
            bl_valid, bl_errors, bl_warnings = tree.validate_branch_lengths()
            if not bl_valid:
                negative_found = any(
                    "Negative" in e or "negative" in e.lower() for e in bl_errors
                )
                if negative_found:
                    raise NegativeBranchLengthError(
                        f"树包含负分支长度: {[e for e in bl_errors if 'Negative' in e or 'negative' in e.lower()]}",
                        suggestion="请检查树文件中分支长度是否正确，或使用 --skip-length-check 跳过此检查。",
                    )
                else:
                    for err in bl_errors:
                        logger.error(err)
            for w in bl_warnings:
                logger.warning(w)
        else:
            logger.info("跳过分支长度检查 (--skip-length-check)")

        # 大规模数据资源提示
        if tree.num_tips > 10000:
            logger.warning(
                f"检测到大型树 ({tree.num_tips} 个末端节点)。"
                "建议：增加线程数 (--threads)、确保充足内存 (>=16GB)、"
                "考虑使用 --low-memory 模式。"
            )

        # DEBUG 级别内存监控
        _log_memory_usage(logger, "树加载完成后")

        # 序列信息
        logger.step(3, "检测序列信息")

        if shutdown_handler.is_shutting_down:
            return EXIT_INTERRUPTED

        from ..infrastructure.alignment_metadata import AlignmentMetadataExtractor

        metadata_extractor = AlignmentMetadataExtractor()
        alignment_metadata = metadata_extractor.extract(parsed_args.sequence)
        logger.info(f"序列格式: {alignment_metadata.format}")

        # §9.2 - PHYLIP格式报ERROR (除非 --ignore-malformed)
        if alignment_metadata.format == "phylip":
            if (
                hasattr(parsed_args, "ignore_malformed")
                and parsed_args.ignore_malformed
            ):
                logger.warning(
                    "PHYLIP格式不被支持，但 --ignore-malformed 已启用，尝试继续。"
                    "建议将序列文件转换为FASTA格式。"
                )
            else:
                logger.error(
                    "PHYLIP格式不被支持。请将序列文件转换为FASTA格式。"
                    "可使用 seqmagick 或 bioawk 工具进行格式转换。"
                )
                return EXIT_DATA_ERROR

        logger.info(f"序列数量: {alignment_metadata.num_sequences}")
        logger.info(f"序列长度: {alignment_metadata.sequence_length} bp")

        seq_type = parsed_args.mol_type or metadata_extractor.validate_sequence_type(
            parsed_args.sequence
        )
        logger.info(f"序列类型: {seq_type}")

        # §9.3 - --no-cross-check: 关闭树与序列末端标签交叉验证
        if not (hasattr(parsed_args, "no_cross_check") and parsed_args.no_cross_check):
            tree_tips = set(tree.tip_names)
            from Bio import SeqIO

            seq_tips = set()
            with open(parsed_args.sequence) as sf:
                for rec in SeqIO.parse(sf, "fasta"):
                    seq_tips.add(rec.id)
            only_in_tree = tree_tips - seq_tips
            only_in_seq = seq_tips - tree_tips
            if only_in_tree or only_in_seq:
                logger.error("树与序列末端标签不一致")
                if only_in_tree:
                    logger.error(f"  仅在树中: {sorted(only_in_tree)[:10]}")
                if only_in_seq:
                    logger.error(f"  仅在序列中: {sorted(only_in_seq)[:10]}")
                return EXIT_DATA_ERROR
            logger.info("树与序列末端标签交叉验证通过")
        else:
            logger.info("交叉验证已通过 --no-cross-check 关闭")

        _log_memory_usage(logger, "序列信息检测后")

        # 加载校准点
        logger.step(4, "解析校准点")

        calibration_loader = CalibrationLoader()
        calibrations = calibration_loader.load(parsed_args.calibrations)
        logger.info(f"加载了 {len(calibrations)} 个校准点")

        # 自动校准：从命令行参数生成校准点并合并
        if parsed_args.auto_calibrate:
            logger.info(f"自动校准参数: {parsed_args.auto_calibrate}")
            auto_cal_strings = parse_auto_calibrate_args(parsed_args.auto_calibrate)
            auto_calibrator = AutoCalibrator(tree)
            calibrated = auto_calibrator.auto_calibrate(auto_cal_strings)
            calibrations.extend(calibrated.calibrations)
            logger.info(f"自动校准后共有 {len(calibrations)} 个校准点")

        if not calibrations:
            logger.error("未加载到任何有效校准点")
            logger.info(
                "请检查校准文件：确保使用支持的约束类型 "
                "(fixed/uniform/soft_lower/maximum/soft_bounds/gamma/skew_normal/skew_t) "
                "以及正确的字段名 (age/min/max/alpha/beta/location/scale/shape/df)"
            )
            return EXIT_DATA_ERROR

        # §4.3 - 分类学表与树标签一致性检查
        if hasattr(parsed_args, "taxonomy_file") and parsed_args.taxonomy_file:
            from ..services.taxonomy_parser import TaxonomyParser

            taxonomy_kwargs = {}
            if (
                hasattr(parsed_args, "taxonomy_delimiter_mode")
                and parsed_args.taxonomy_delimiter_mode
            ):
                taxonomy_kwargs["delimiter_mode"] = parsed_args.taxonomy_delimiter_mode
            if (
                hasattr(parsed_args, "taxonomy_source_priority")
                and parsed_args.taxonomy_source_priority
            ):
                taxonomy_kwargs["source_priority"] = (
                    parsed_args.taxonomy_source_priority
                )
            if hasattr(parsed_args, "taxonomy_levels") and parsed_args.taxonomy_levels:
                levels = {}
                for item in parsed_args.taxonomy_levels:
                    if ":" in item:
                        prefix, rank = item.split(":", 1)
                        levels[prefix] = rank
                    else:
                        levels[item] = item
                if levels:
                    taxonomy_kwargs["custom_levels"] = levels

            taxonomy_parser = TaxonomyParser(**taxonomy_kwargs)

            load_kwargs = {}
            if hasattr(parsed_args, "table_sep") and parsed_args.table_sep:
                load_kwargs["table_sep"] = parsed_args.table_sep
            if (
                hasattr(parsed_args, "taxonomy_table_sep")
                and parsed_args.taxonomy_table_sep
            ):
                load_kwargs["taxonomy_sep"] = parsed_args.taxonomy_table_sep

            # 审阅项 B-26 [已证]：旧实现丢弃了 ``load_from_file()`` 的返回值，
            # 并且无条件解包 ``check_label_consistency()``。分类学表整份加载
            # 失败时（表头缺名称列、全是解析不出的行……）一致性检查会把**全部**
            # 叶节点报成"已匹配"、两个未匹配清单为空，于是这项质量检查在表
            # 完全失效时呈现为 100% 通过。现在：加载失败即数据错误退出；
            # 检查未执行（返回 None）同样退出；零交集也视为错误。
            try:
                loaded_entries = taxonomy_parser.load_from_file(
                    parsed_args.taxonomy_file, **load_kwargs
                )
            except (FileNotFoundError, ValueError) as e:
                logger.error(f"分类学表加载失败: {parsed_args.taxonomy_file}: {e}")
                logger.info(
                    "多列表格必须含名称列 (name/taxon/organism/label)，"
                    "两列表格第二列须是完整分类串；"
                    "整表加载失败时不得继续，否则类群↔化石锚定关系无从核对"
                )
                return EXIT_DATA_ERROR

            logger.info(
                f"已加载 {loaded_entries} 条分类学条目: {parsed_args.taxonomy_file}"
            )

            consistency = taxonomy_parser.check_label_consistency(tree.tip_names)
            if consistency is None:
                logger.error(
                    "分类学表↔树标签一致性检查未能执行："
                    f"{parsed_args.taxonomy_file} 中没有可用于比对的条目"
                )
                logger.info(
                    "未执行 ≠ 通过。请修正表格（表头含名称列与分类级别列，"
                    "或第二列为完整分类串）后重跑"
                )
                return EXIT_DATA_ERROR

            matched, unmatched_tree, unmatched_table = consistency
            if not matched:
                logger.error(
                    f"分类学表与树标签零匹配：{len(tree.tip_names)} 个末端节点中"
                    f"没有任何一个出现在 {loaded_entries} 条分类学条目里"
                    f"（表中未匹配条目: {unmatched_table[:5]}）"
                )
                logger.info(
                    "零匹配通常意味着表格与这棵树并非同一批样本，"
                    "或命名格式（分隔符/等级前缀）不一致；"
                    "可用 --table-sep / --taxonomy-table-sep / "
                    "--taxonomy-delimiter-mode 指定后重试"
                )
                return EXIT_DATA_ERROR

            logger.info(
                f"分类学表↔树标签一致性检查: {len(matched)}/"
                f"{len(tree.tip_names)} 个末端节点已匹配"
            )
            if unmatched_tree:
                logger.warning(
                    f"{len(unmatched_tree)} tree tips not found in taxonomy table: "
                    f"{unmatched_tree[:5]}..."
                )
            if unmatched_table:
                logger.warning(
                    f"{len(unmatched_table)} taxonomy entries not found in tree: "
                    f"{unmatched_table[:5]}..."
                )

        # 应用方法专有参数到配置
        config = Configuration()
        if parsed_args.config:
            try:
                config.load_from_file(parsed_args.config)
            except Exception as e:
                logger.error(f"配置文件解析失败: {e}")
                return EXIT_USAGE_ERROR

        # §8.1 - CLI参数优先于Config
        _apply_method_args(config, parsed_args.method, parsed_args.method_args, logger)

        # 构建 PipelineConfig
        # 注意：burnin/nsample/timeout 若通过 --method-args 以 mcmctree.burnin 等形式
        # 设置，需显式透传到 PipelineConfig，否则 pipeline.__init__ 会用 PipelineConfig
        # 默认值(20000/50000)覆盖 tool_config 中用户已通过 method-args 设定的值(F15)。
        # 未设置时不传(而非传 None)，让 PipelineConfig 保留自身默认值。
        _pipeline_kwargs: dict = {}
        if parsed_args.method_args:
            _ma_burnin = parsed_args.method_args.get("burnin")
            _ma_nsample = parsed_args.method_args.get("nsample")
            _ma_timeout = parsed_args.method_args.get("timeout")
            if _ma_burnin is not None:
                _pipeline_kwargs["burnin"] = int(_ma_burnin)
            if _ma_nsample is not None:
                _pipeline_kwargs["nsample"] = int(_ma_nsample)
            if _ma_timeout is not None:
                _pipeline_kwargs["timeout"] = int(_ma_timeout)
        # 构建 taxonomy 相关参数
        _taxonomy_kwargs: dict = {}
        if (
            hasattr(parsed_args, "taxonomy_delimiter_mode")
            and parsed_args.taxonomy_delimiter_mode
        ):
            _taxonomy_kwargs["taxonomy_delimiter_mode"] = (
                parsed_args.taxonomy_delimiter_mode
            )
        if (
            hasattr(parsed_args, "taxonomy_source_priority")
            and parsed_args.taxonomy_source_priority
        ):
            _taxonomy_kwargs["taxonomy_source_priority"] = (
                parsed_args.taxonomy_source_priority
            )
        if hasattr(parsed_args, "taxonomy_levels") and parsed_args.taxonomy_levels:
            _levels = {}
            for item in parsed_args.taxonomy_levels:
                if ":" in item:
                    _prefix, _rank = item.split(":", 1)
                    _levels[_prefix] = _rank
                else:
                    _levels[item] = item
            if _levels:
                _taxonomy_kwargs["taxonomy_levels"] = _levels
        if hasattr(parsed_args, "table_sep") and parsed_args.table_sep:
            _taxonomy_kwargs["taxonomy_table_sep"] = parsed_args.table_sep
        if (
            hasattr(parsed_args, "taxonomy_table_sep")
            and parsed_args.taxonomy_table_sep
        ):
            _taxonomy_kwargs["taxonomy_sep"] = parsed_args.taxonomy_table_sep

        pipeline_config = PipelineConfig(
            methods=[parsed_args.method],
            output_dir=output_dir,
            threads=parsed_args.threads,
            seed=parsed_args.seed,
            dry_run=parsed_args.dry_run,
            software_paths=config.software_paths,
            taxonomy_file=parsed_args.taxonomy_file,
            reroot_strategy=parsed_args.reroot,
            outgroup_name=parsed_args.outgroup,
            paml_version=parsed_args.method_args.get("paml_version"),
            **_pipeline_kwargs,
            **_taxonomy_kwargs,
        )

        # §11.1 - --resume: 从检查点恢复
        if hasattr(parsed_args, "resume") and parsed_args.resume:
            pipeline_config.resume_mode = True
            logger.info("Resume mode: will skip completed checkpoint steps")

        # §16.1 - --low-memory: 传入Pipeline
        if hasattr(parsed_args, "low_memory") and parsed_args.low_memory:
            pipeline_config.low_memory = True
            logger.info("Low-memory mode enabled")

        # --strip-annotations: 输出时丢弃树注释
        pipeline_config.strip_annotations = parsed_args.strip_annotations
        if pipeline_config.strip_annotations:
            logger.info("Strip annotations mode: tree comments will be discarded")

        logger.step(5, "执行定年分析")

        if parsed_args.threads > 1:
            # ParallelPipeline 是 DatingPipeline 的子类；先窄后宽地复用同一个名字，
            # 必须把变量标成公共类型，否则第二个赋值会被判为"子类→父类"的错。
            pipeline: DatingPipeline
            pipeline = ParallelPipeline(
                pipeline_config, tool_config=config.get_config()
            )
        else:
            pipeline = DatingPipeline(pipeline_config, tool_config=config.get_config())

        # 将 pipeline 的检查点管理器注入 shutdown handler，以便中断时保存进度
        shutdown_handler._state.checkpoint_manager = pipeline._checkpoint_manager

        results = pipeline.run(tree, parsed_args.sequence, calibrations)

        # §14 - 创建运行标记
        output_marker.write_text("completed")

        # 汇总
        logger.section("定年分析完成")
        elapsed = time.time() - start_time
        if parsed_args.dry_run:
            # dry-run 不执行真实分析，results 为空属正常，避免"完成方法数: 0"
            # 让用户误以为什么都没干。
            logger.info("Dry-run 模式：配置与依赖校验通过，未执行真实分析")
        else:
            logger.info(f"完成方法数: {len(results)}")
        logger.timing("总耗时", elapsed)

        if parsed_args.dry_run:
            return EXIT_SUCCESS
        return EXIT_SUCCESS if results else EXIT_RUNTIME_ERROR

    except KeyboardInterrupt:
        logger.warning("\n用户中断")
        return EXIT_INTERRUPTED

    except SystemExit as e:
        if e.code == 2:
            return EXIT_USAGE_ERROR
        raise

    except (
        TreeValidationError,
        AlignmentValidationError,
        ConsistencyError,
        CalibrationError,
    ) as e:
        logger.error(f"数据格式错误: {e}")
        if hasattr(e, "suggestion") and e.suggestion:
            logger.info(f"建议: {e.suggestion}")
        return EXIT_DATA_ERROR

    except (ConfigurationError, UnknownMethodError, ValueError) as e:
        logger.error(f"配置或参数错误: {e}")
        if hasattr(e, "suggestion") and e.suggestion:
            logger.info(f"建议: {e.suggestion}")
        return EXIT_DATA_ERROR

    except Exception as e:
        logger.error(f"运行错误: {e}")
        import traceback

        logger.debug_file_only(traceback.format_exc())
        return EXIT_RUNTIME_ERROR

    finally:
        shutdown_handler.unregister()
        logger.stop()


def _apply_method_args(
    config: Configuration,
    method: str,
    method_args: Dict[str, str],
    logger: "PhyloDaterLogger",
) -> None:
    """将 --method-args 键值对应用到配置"""
    if not method_args:
        return

    # 通用参数映射
    common_mappings = {
        "paml_path": ("software_paths", "paml_path"),
        "mcmctree_bin": ("software_paths", "mcmctree_bin"),
        "iqtree_bin": ("software_paths", "iqtree_bin"),
    }

    for key, value in method_args.items():
        # 检查是否是通用参数
        if key in common_mappings:
            attr_path = common_mappings[key]
            config.set_software_path(attr_path[1], value)
            logger.info(f"设置 {key} = {value}")
            continue

        # 作为方法专有参数传递
        config.set_cli_override(method, key, value)
        logger.info(f"设置 {method}.{key} = {value}")


def _run_fold(parsed_args: argparse.Namespace, logger: "PhyloDaterLogger") -> int:
    """执行树折叠子命令。

    该子命令尚未实现。**先**打出"未实现"警告再回显参数（审阅项 C-36）：旧实现
    把解析结果一路回显到最后才 raise，用户读到的前几行与已实现的子命令毫无区别。
    随后显式抛出 NotImplementedError，让用户拿到明确错误与非零退出码，而不是
    返回成功（假阳性）或写出一个"看起来折叠过、其实没动"的树文件。
    调用方（main）捕获后返回 EXIT_RUNTIME_ERROR。
    """
    logger.section("PhyloDater 树折叠")
    logger.warning(
        f"fold 子命令【{EXPERIMENTAL_NOT_IMPLEMENTED}】：以下参数仅为接口预留，"
        "本次调用不会产出任何折叠树。"
    )
    logger.info(f"输入树: {parsed_args.tree}")
    logger.info(f"折叠级别: {parsed_args.fold_rank}")
    raise NotImplementedError(
        "fold 子命令尚未实现（--tree/--output/--fold-rank 均为接口预留参数）"
    )


# ==================== 主入口 ====================


def main(args: Optional[List[str]] = None) -> int:
    """主入口函数"""

    # ------------------------------------------------------------------
    # 预解析通用参数（--log-file/--log-level 等），使其既可放在子命令前，
    # 也可放在子命令后。argparse 的 parent parser 机制会把子命令前的参数
    # 仅解析到主解析器，而子命令的 namespace 中却丢失该值，因此这里用
    # 一个独立的 common_parser 扫描整组参数，再与完整解析结果合并。
    # ------------------------------------------------------------------
    common_parser = argparse.ArgumentParser(add_help=False)
    _add_common_arguments(common_parser)
    common_ns, remaining = common_parser.parse_known_args(args)

    parser = create_parser()

    try:
        parsed_args = parser.parse_args(remaining)
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else EXIT_USAGE_ERROR

    # 将 common 参数合并到最终 namespace 中（支持前后两种位置）
    _COMMON_ATTRS = [
        "log_file",
        "log_level",
        "taxonomy_file",
        "taxonomy_levels",
        "taxonomy_delimiter_mode",
        "taxonomy_source_priority",
        "table_sep",
        "taxonomy_table_sep",
        "ignore_malformed",
    ]
    for attr in _COMMON_ATTRS:
        common_val = getattr(common_ns, attr, None)
        if common_val is not None:
            setattr(parsed_args, attr, common_val)

    # §17 - 钳制 --threads 到 [1, cpu_count]（负值已在前端拒绝）
    if hasattr(parsed_args, "threads"):
        parsed_args.threads = _clamp_threads(parsed_args.threads)

    show_banner()

    # 检查是否指定了子命令
    if not parsed_args.subcommand:
        parser.print_help()
        return EXIT_USAGE_ERROR

    # 子命令路由
    if parsed_args.subcommand in ("check", "self-test"):
        # 自检模式：不需要输出目录
        import tempfile

        output_dir = Path(tempfile.mkdtemp(prefix="phylodater_check_"))
        # 显式设置临时目录权限为 0o700，确保仅当前用户可访问
        try:
            os.chmod(str(output_dir), 0o700)
        except Exception:
            pass
        logger = _setup_logger(parsed_args, output_dir)

        from ..infrastructure.signal_handler import get_shutdown_handler

        shutdown_handler = get_shutdown_handler()
        shutdown_handler.register()

        try:
            return _run_self_test(logger)
        finally:
            shutdown_handler.unregister()
            logger.stop()
            import shutil

            shutil.rmtree(output_dir, ignore_errors=True)

    elif parsed_args.subcommand == "dating":
        output_dir = Path(parsed_args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger = _setup_logger(parsed_args, output_dir)
        return _run_dating(parsed_args, logger)

    elif parsed_args.subcommand == "fold":
        output_dir = Path(parsed_args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger = _setup_logger(parsed_args, output_dir)
        try:
            return _run_fold(parsed_args, logger)
        except NotImplementedError as e:
            # fold 子命令尚未实现：显式给出错误而非返回成功（假阳性）
            logger.error(f"功能未实现: {e}")
            logger.error(
                "该能力已在 --help 中标注为实验性；本版本不提供树折叠功能，"
                "输出目录中不会有任何结果文件。"
            )
            return EXIT_RUNTIME_ERROR
        finally:
            logger.stop()

    else:
        print(f"未知子命令: {parsed_args.subcommand}")
        parser.print_help()
        return EXIT_USAGE_ERROR


def _validate_input_files(
    parsed_args: argparse.Namespace, logger: "PhyloDaterLogger"
) -> bool:
    """验证输入文件是否有效"""
    validator = InputValidator()
    ignore_malformed = getattr(parsed_args, "ignore_malformed", False)
    deep_validator = DeepValidator(ignore_malformed=ignore_malformed)
    all_valid = True

    # 验证树文件
    logger.info("验证树文件...")
    tree_result = validator.validate_tree_file(parsed_args.tree)
    if tree_result.is_valid:
        logger.success(f"树文件格式有效: {tree_result.format_detected}")
        for info in tree_result.info:
            logger.info(info)
        for warning in tree_result.warnings:
            logger.warning(warning)
    else:
        for error in tree_result.errors:
            logger.error(error)
        all_valid = False

    # 深度验证树文件
    logger.info("深度验证树文件结构...")
    tree_detail = deep_validator.validate_tree_deep(
        parsed_args.tree, ignore_malformed=ignore_malformed
    )
    if tree_detail.is_valid:
        logger.success("树文件深度验证通过")
        if tree_detail.tree_count > 1:
            logger.warning(f"树文件包含 {tree_detail.tree_count} 棵树")
            for summary in tree_detail.tree_summaries:
                tree_name = summary.get("name", f"树 {summary['tree_num']}")
                logger.info(
                    f"  {tree_name}: {summary['terminals']} 个末端节点, {summary['internal']} 个内部节点"
                )
        for warning in tree_detail.warnings:
            logger.warning(f"  {warning}")
    else:
        logger.error("树文件深度验证失败:")
        for error in tree_detail.errors:
            logger.error(f"  {error}")
        for warning in tree_detail.warnings:
            logger.warning(f"  {warning}")
        all_valid = False

    # 验证序列文件
    logger.info("验证序列比对文件...")
    alignment_result = validator.validate_alignment_file(parsed_args.sequence)
    if alignment_result.is_valid:
        logger.success(f"序列比对文件格式有效: {alignment_result.format_detected}")
        for info in alignment_result.info:
            logger.info(info)
        for warning in alignment_result.warnings:
            logger.warning(warning)
    else:
        for error in alignment_result.errors:
            logger.error(error)
        all_valid = False

    # 深度验证序列文件
    logger.info("深度验证序列文件...")
    seq_detail = deep_validator.validate_alignment_deep(
        parsed_args.sequence,
        expected_alphabet=None,
        require_aligned=not parsed_args.skip_length_check,
        ignore_malformed=ignore_malformed,
    )
    if seq_detail.is_valid:
        logger.success("序列文件深度验证通过")
        logger.info(f"  序列数量: {seq_detail.num_sequences}")
        if seq_detail.sequence_length:
            logger.info(f"  序列长度: {seq_detail.sequence_length} bp")
        if seq_detail.alphabet_detected:
            logger.info(f"  字母表类型: {seq_detail.alphabet_detected.upper()}")
        if not seq_detail.is_aligned:
            logger.warning("  序列长度不一致（未对齐）")
        for warning in seq_detail.warnings:
            logger.warning(f"  {warning}")
    else:
        logger.error("序列文件深度验证失败:")
        for error in seq_detail.errors:
            logger.error(f"  {error}")
        for warning in seq_detail.warnings:
            logger.warning(f"  {warning}")
        all_valid = False

    # 验证校准配置文件
    logger.info("验证校准配置文件...")
    calibration_result = validator.validate_calibration_file(parsed_args.calibrations)
    if calibration_result.is_valid:
        logger.success(f"校准配置文件格式有效: {calibration_result.format_detected}")
        for warning in calibration_result.warnings:
            logger.warning(warning)
    else:
        for error in calibration_result.errors:
            logger.error(error)
        all_valid = False

    # 验证分类学文件（如果提供）
    if hasattr(parsed_args, "taxonomy_file") and parsed_args.taxonomy_file:
        logger.info("验证分类学信息文件...")
        taxonomy_result = validator.validate_taxonomy_file(parsed_args.taxonomy_file)
        if taxonomy_result.is_valid:
            logger.success(f"分类学信息文件格式有效: {taxonomy_result.format_detected}")
            for warning in taxonomy_result.warnings:
                logger.warning(warning)
        else:
            for error in taxonomy_result.errors:
                logger.error(error)
            all_valid = False

    return all_valid


def _log_memory_usage(logger: "PhyloDaterLogger", context: str = "") -> None:
    """DEBUG 级别记录内存使用（psutil 可选）"""
    try:
        import psutil

        process = psutil.Process()
        mem_mb = process.memory_info().rss / (1024 * 1024)
        sys_mem = psutil.virtual_memory()
        logger.debug(
            f"内存使用 {context}: 进程 {mem_mb:.1f} MB, "
            f"系统可用 {sys_mem.available / (1024**3):.1f} GB / "
            f"总计 {sys_mem.total / (1024**3):.1f} GB"
        )
    except ImportError:
        pass  # psutil 不可用时忽略


if __name__ == "__main__":
    sys.exit(main())
