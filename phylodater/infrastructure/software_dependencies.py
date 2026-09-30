"""
SoftwareDependencyManager - 软件依赖管理器

提供各定年软件的安装建议，支持跨平台
"""

import platform
import shutil
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional


@dataclass
class SoftwareInfo:
    """软件信息"""

    name: str
    executable: str
    conda_package: str
    conda_channel: str
    brew_package: Optional[str] = None
    apt_package: Optional[str] = None
    description: str = ""


class ProbeStatus(str, Enum):
    """可执行文件探测结果（审阅项 C-57）。

    三种状态必须分开：把"探测过程本身失败"与"确认未安装"塌缩成同一个
    ``False``，用户看到的效果是"照提示重装一遍还是缺失"，而真实原因可能在
    别处（PATH 目录不可读、权限、名字拼错）。本枚举让上层至少能问出
    "到底是没装，还是我根本没查成"。
    """

    AVAILABLE = "available"  # 正向定位到可执行文件
    MISSING = "missing"  # 探测成功，但 PATH 上确实没有
    PROBE_FAILED = "probe_failed"  # 探测过程出错，可用性**未知**


@dataclass
class DependencyProbe:
    """单个软件的探测结果。"""

    software: str
    status: ProbeStatus
    path: Optional[str] = None
    error: Optional[str] = None

    @property
    def available(self) -> bool:
        """仅在**正向定位到**可执行文件时为 True；探测失败按未知处理。"""
        return self.status is ProbeStatus.AVAILABLE

    def __str__(self) -> str:
        if self.status is ProbeStatus.AVAILABLE:
            return f"{self.software}: available ({self.path})"
        if self.status is ProbeStatus.MISSING:
            return f"{self.software}: not found in PATH"
        return f"{self.software}: probe failed ({self.error})"


SOFTWARE_CATALOG: Dict[str, SoftwareInfo] = {
    "mcmctree": SoftwareInfo(
        name="MCMCTree",
        executable="mcmctree",
        conda_package="paml",
        conda_channel="bioconda",
        apt_package="paml",
        description="Bayesian molecular dating tool from PAML package",
    ),
    "treepl": SoftwareInfo(
        name="treePL",
        executable="treePL",
        conda_package="treepl",
        conda_channel="bioconda",
        apt_package=None,
        description="Maximum likelihood divergence time estimation",
    ),
    "pathd8": SoftwareInfo(
        name="PATHd8",
        executable="PATHd8",
        conda_package="pathd8",
        conda_channel="bioconda",
        apt_package=None,
        description="Phylogenetic dating tool for non-clock trees",
    ),
    "iqtree2": SoftwareInfo(
        name="IQ-TREE2",
        executable="iqtree2",
        conda_package="iqtree2",
        conda_channel="bioconda",
        brew_package="iqtree",
        apt_package="iqtree2",
        description="Fast and effective IQ-TREE2 phylogenetic inference",
    ),
    "baseml": SoftwareInfo(
        name="BASEML",
        executable="baseml",
        conda_package="paml",
        conda_channel="bioconda",
        apt_package="paml",
        description="BASEML tool from PAML package",
    ),
    "codeml": SoftwareInfo(
        name="CODEML",
        executable="codeml",
        conda_package="paml",
        conda_channel="bioconda",
        apt_package="paml",
        description="CODEML tool from PAML package",
    ),
    "pyr8s": SoftwareInfo(
        name="pyr8s",
        executable="pyr8s",
        conda_package="pyr8s",
        conda_channel="bioconda",
        apt_package=None,
        description="Python wrapper for r8s dating tool",
    ),
    "mdcat": SoftwareInfo(
        name="MD-Cat",
        # 可执行/包名核对（审阅项 C-56 附带项）：上游仓库
        # `PhyloDater-参考软件/原始代码-MD-Cat` 的入口脚本是 `md_cat.py`，
        # setup.py 以 `scripts=['md_cat.py', 'simulate.py']` 安装（pip 装完后
        # 命令名保留 .py 后缀），发行名是 `MD-Cat`（见 MD_Cat.egg-info/PKG-INFO：
        # `Name: MD-Cat`、`Version: 1.0.1`）。也就是说上游**没有** `mdcat` 这个
        # 可执行名。anaconda.org 检索 `mdcat` 只命中 conda-forge 的
        # "Show markdown documents on text terminals"（终端 Markdown 查看器，与本
        # 软件无关），未见于 bioconda；而 conda 的常规配置里 bioconda 依赖
        # conda-forge 通道，于是本目录生成的 `conda install -c bioconda mdcat`
        # 很可能装上那个错的东西（此推断待人工复核，见交付说明）。
        # 适配器实际探测/调用的路径另有其来源：配置项 `software_paths.mdcat_bin`
        # （`adapters/mdcat_method.py:80-83`），其代码默认值是 `mdcat`
        # （`infrastructure/configuration.py:374`），config.example.yaml 里却写
        # `MD-Cat` —— 三处命名互不相同。把它们对齐需要改 adapters/、
        # infrastructure/configuration.py 与文档，属跨文件改动，本轮不动。
        executable="mdcat",
        conda_package="mdcat",
        conda_channel="bioconda",
        apt_package=None,
        # 原名解释有误："Cat" 指速率**类别**（categorical rate distribution），
        # 与 Catmull-Rom/Friedman 样条定年无关。判据见上游 README：
        # "MD-Cat relaxes the molecular clock assumption by approximating the rate
        #  distribution by a categorical distribution. The only required parameter
        #  for this model is the number of rate categories, which is default to 50."
        # 以及 PKG-INFO 的 Summary: "Phylogenetic dating under a flexible
        # categorical model using Expectation-Maximization"。
        description=(
            "Molecular Dating with rate Categories " "(categorical rate distribution)"
        ),
    ),
    "wlogdate": SoftwareInfo(
        name="wLogDate",
        executable="launch_wLogDate.py",
        conda_package="wlogdate",
        conda_channel="bioconda",
        apt_package=None,
        description="Weighted log dating for molecular phylogenies",
    ),
}


class SoftwareDependencyManager:
    """
    软件依赖管理器

    提供缺失软件的具体安装建议
    """

    @staticmethod
    def get_install_command(executable: str) -> str:
        """
        获取指定软件的安装命令

        Args:
            executable: 可执行文件名（如 'mcmctree', 'treePL'）

        Returns:
            针对当前操作系统的具体安装命令
        """
        system = platform.system()
        executable_lower = executable.lower().replace("-", "").replace("_", "")

        info = None
        for key, software_info in SOFTWARE_CATALOG.items():
            if (
                key.replace("_", "") == executable_lower
                or software_info.executable.lower().replace("-", "").replace("_", "")
                == executable_lower
            ):
                info = software_info
                break

        if info is None:
            return f"# Software '{executable}' not found in catalog. Please install manually."

        if system == "Darwin":
            return SoftwareDependencyManager._get_macos_install(info)
        elif system == "Linux":
            return SoftwareDependencyManager._get_linux_install(info)
        else:
            return f"# Unsupported platform: {system}. Please install {info.name} manually."

    @staticmethod
    def _get_macos_install(info: SoftwareInfo) -> str:
        """生成 macOS 安装命令"""
        commands = []
        commands.append(f"# Install {info.name} on macOS")
        commands.append("")
        commands.append("# Option 1: Using conda (recommended)")
        if info.conda_channel == "bioconda":
            commands.append(f"conda install -c bioconda {info.conda_package}")
        else:
            commands.append(
                f"conda install -c {info.conda_channel} {info.conda_package}"
            )
        commands.append("")
        if info.brew_package:
            commands.append("# Option 2: Using Homebrew")
            commands.append(f"brew install {info.brew_package}")
        commands.append("")
        commands.append("# Option 3: Using MacPorts")
        if info.apt_package:
            commands.append(f"# sudo port install {info.apt_package}")
        return "\n".join(commands)

    @staticmethod
    def _get_linux_install(info: SoftwareInfo) -> str:
        """生成 Linux 安装命令"""
        commands = []
        commands.append(f"# Install {info.name} on Linux")
        commands.append("")
        commands.append("# Option 1: Using conda (recommended)")
        if info.conda_channel == "bioconda":
            commands.append(f"conda install -c bioconda {info.conda_package}")
        else:
            commands.append(
                f"conda install -c {info.conda_channel} {info.conda_package}"
            )
        commands.append("")
        if info.apt_package:
            commands.append("# Option 2: Using apt-get")
            commands.append(f"sudo apt-get install {info.apt_package}")
        commands.append("")
        commands.append("# Option 3: Build from source")
        commands.append(f"# See {info.name} documentation for build instructions")
        return "\n".join(commands)

    @staticmethod
    def probe_dependencies(required_software: List[str]) -> Dict[str, DependencyProbe]:
        """逐个探测可执行文件，返回三态结果（审阅项 C-57）。

        用 :func:`shutil.which` 而非 ``subprocess.run(["which", …])``：
        ``which`` 是 POSIX 外部命令，Windows 默认没有它，旧实现在那里会抛
        ``FileNotFoundError`` 并被 ``except Exception`` 吞掉 → 每个依赖都被
        报成"未安装"。``shutil.which`` 是纯 Python 实现、跨平台（Windows 上
        还会处理 ``PATHEXT``/.exe），并且不起子进程。本仓库其余各处
        （``process_runner.py``、``container.py``、``treepl_method.py``、
        ``mcmctree_method.py``）早已用它，这里是回到统一正解。

        Args:
            required_software: 待探测的可执行文件名/软件名列表

        Returns:
            {软件名: DependencyProbe}，其中 status 区分
            ``AVAILABLE`` / ``MISSING`` / ``PROBE_FAILED``
        """
        probes: Dict[str, DependencyProbe] = {}
        for software in required_software:
            try:
                path = shutil.which(software)
            except OSError as e:
                # shutil.which 只在 PATH 含无法访问的目录等极端情况下抛 OSError；
                # 这属于"探测没做成"，不等于"软件没装"。
                probes[software] = DependencyProbe(
                    software=software,
                    status=ProbeStatus.PROBE_FAILED,
                    error=f"{type(e).__name__}: {e}",
                )
                continue
            if path:
                probes[software] = DependencyProbe(
                    software=software, status=ProbeStatus.AVAILABLE, path=path
                )
            else:
                probes[software] = DependencyProbe(
                    software=software, status=ProbeStatus.MISSING
                )
        return probes

    @staticmethod
    def check_all_dependencies(required_software: list) -> Dict[str, bool]:
        """
        检查所有依赖是否可用（布尔便捷视图）

        Args:
            required_software: 所需软件列表

        Returns:
            {软件名: 是否可用}

        Note:
            ``True`` 只在**正向定位到**可执行文件时给出。探测失败（例如 PATH
            目录不可读）在本视图里保守地记为 ``False``（fail-closed：绝不把
            "未知"说成"可用"），但它与"确实未安装"是两回事，会额外记一条
            WARNING；需要区分二者的调用方应改用
            :meth:`probe_dependencies` 并读取 ``status``。
        """
        probes = SoftwareDependencyManager.probe_dependencies(required_software)

        probe_failures = [
            probe
            for probe in probes.values()
            if probe.status is ProbeStatus.PROBE_FAILED
        ]
        if probe_failures:
            try:
                from .logging import get_logger

                get_logger().warning(
                    "Dependency probe failed (availability UNKNOWN, not 'missing') "
                    "for: "
                    + ", ".join(f"{p.software} [{p.error}]" for p in probe_failures)
                )
            except Exception:  # pragma: no cover - 日志不可用时不影响返回值
                pass

        return {software: probe.available for software, probe in probes.items()}
