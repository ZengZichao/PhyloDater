"""
Container Environment - 容器化环境兼容

提供容器化环境检测和路径管理功能，确保管线在以下环境中正确运行：
- Docker
- Singularity
- HPC集群 (PBS, SLURM)
- 本地开发环境
"""

import os
import tempfile
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Optional


class ContainerType(Enum):
    """容器类型"""

    NONE = "none"
    DOCKER = "docker"
    SINGULARITY = "singularity"
    PBS = "pbs"
    SLURM = "slurm"


@dataclass
class ContainerEnvironment:
    """
    容器环境信息

    检测并记录当前运行环境类型和可用资源
    """

    container_type: ContainerType
    has_write_permission: bool
    temp_dir: Path
    work_dir: Path
    cache_dir: Optional[Path] = None
    #: ``dataclass`` 的默认值不能是 ``None`` 却标成 ``Dict``：旧写法靠
    #: ``__post_init__`` 里补 ``{}`` 才能跑，但静态视角下那个字段永远是
    #: ``Dict[str, Any]``，给它的赋值就都成了类型错。``default_factory`` 把
    #: 真实行为（默认空字典、允许调用方显式传 None 补平）写进类型里。
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.metadata is None:
            self.metadata = {}


def detect_container_type() -> ContainerType:
    """
    检测当前运行的容器类型

    Returns:
        ContainerType 枚举值
    """
    if os.environ.get("SINGULARITY_CONTAINER"):
        return ContainerType.SINGULARITY
    elif os.environ.get("DOCKER"):
        return ContainerType.DOCKER
    elif os.environ.get("PBS_ENVIRONMENT"):
        return ContainerType.PBS
    elif os.environ.get("SLURM_JOB_ID"):
        return ContainerType.SLURM
    else:
        return ContainerType.NONE


def get_writable_directory(preferred: Optional[Path] = None) -> Path:
    """
    获取具有写权限的目录

    在容器化环境中，优先使用以下目录：
    - /tmp (Docker 默认)
    - $SINGULARITY_CONTAINER_DIR (Singularity)
    - 用户主目录下的 .phylodater 目录

    Args:
        preferred: 优先使用的目录路径

    Returns:
        可写的目录路径
    """
    if preferred is not None:
        preferred = Path(preferred)
        if os.access(preferred, os.W_OK):
            return preferred
        if not preferred.exists():
            try:
                preferred.mkdir(parents=True, exist_ok=True)
                return preferred
            except PermissionError:
                pass

    container_type = detect_container_type()

    if container_type == ContainerType.SINGULARITY:
        singularity_dir = os.environ.get("SINGULARITY_CONTAINER_DIR")
        if singularity_dir:
            singularity_path = Path(singularity_dir)
            if os.access(singularity_path, os.W_OK):
                return singularity_path

    temp_base = tempfile.gettempdir()
    temp_dir = Path(temp_base)

    if os.access(temp_dir, os.W_OK):
        phylodater_temp = temp_dir / ".phylodater"
        try:
            phylodater_temp.mkdir(exist_ok=True)
            return phylodater_temp
        except PermissionError:
            return temp_dir

    home_dir = Path.home()
    phylodater_home = home_dir / ".phylodater"
    try:
        phylodater_home.mkdir(exist_ok=True)
        return phylodater_home
    except PermissionError:
        return home_dir


def get_container_specific_env() -> Dict[str, str]:
    """
    获取容器特定的环境变量

    Returns:
        容器相关的环境变量字典
    """
    env_vars: Dict[str, str] = {}

    container_type = detect_container_type()

    if container_type == ContainerType.SINGULARITY:
        singularity_env = [
            "SINGULARITY_CONTAINER",
            "SINGULARITY_CONTAINER_DIR",
            "SINGULARITY_NAME",
            "SINGULARITY_APPNAME",
        ]
        for var in singularity_env:
            value = os.environ.get(var)
            if value:
                env_vars[var] = value

    elif container_type == ContainerType.DOCKER:
        docker_env = [
            "DOCKER_CONTAINER",
        ]
        for var in docker_env:
            value = os.environ.get(var)
            if value:
                env_vars[var] = value

    elif container_type == ContainerType.SLURM:
        slurm_env = [
            "SLURM_JOB_ID",
            "SLURM_JOB_NODELIST",
            "SLURM_NTASKS",
            "SLURM_SUBMIT_DIR",
        ]
        for var in slurm_env:
            value = os.environ.get(var)
            if value:
                env_vars[var] = value

    return env_vars


def check_software_in_container(software_name: str) -> Optional[str]:
    """
    检查容器内是否安装了指定的软件

    Args:
        software_name: 软件名称

    Returns:
        软件路径或 None
    """
    import shutil

    path = shutil.which(software_name)
    if path:
        return path

    common_paths = [
        Path("/usr/local/bin") / software_name,
        Path("/usr/bin") / software_name,
        Path("/opt") / software_name / software_name,
    ]

    if detect_container_type() == ContainerType.SINGULARITY:
        container_dir = os.environ.get("SINGULARITY_CONTAINER_DIR", "")
        if container_dir:
            common_paths.insert(0, Path(container_dir) / "bin" / software_name)

    # 旧实现在这里把循环变量也叫 ``path``，于是同一个名字先后承着
    # ``Optional[str]``（shutil.which 的结果）与 ``Path``（候选路径）两种类型。
    for candidate in common_paths:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)

    return None


def create_container_compatible_workdir(
    base_dir: Optional[Path] = None, prefix: str = "phylodater"
) -> Path:
    """
    创建容器兼容的工作目录

    Args:
        base_dir: 基础目录，None 则使用可写目录
        prefix: 目录名前缀

    Returns:
        创建的工作目录路径
    """
    import time
    import uuid

    if base_dir is None:
        base_dir = get_writable_directory()

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    unique_id = uuid.uuid4().hex[:8]
    work_dir = base_dir / f"{prefix}_{timestamp}_{unique_id}"

    try:
        work_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        fallback_dir = get_writable_directory()
        work_dir = fallback_dir / f"{prefix}_{timestamp}_{unique_id}"
        work_dir.mkdir(parents=True, exist_ok=True)

    return work_dir


def get_container_environment_info() -> ContainerEnvironment:
    """
    获取完整的容器环境信息

    Returns:
        ContainerEnvironment 对象
    """
    container_type = detect_container_type()
    work_dir = get_writable_directory()
    temp_dir = Path(tempfile.gettempdir())

    env_info = get_container_specific_env()

    # 只有建不起来的目录才变成 None：先把路径建出来，再决定要不要记住它。
    # （旧写法把同一个变量先当 Path 用、再当 None 存，静态视角下它一路都是
    # Optional，mkdir 就成了 union-attr。）
    cache_dir: Optional[Path]
    cache_path = work_dir / ".cache"
    try:
        cache_path.mkdir(exist_ok=True)
        cache_dir = cache_path
    except PermissionError:
        cache_dir = None

    return ContainerEnvironment(
        container_type=container_type,
        has_write_permission=os.access(work_dir, os.W_OK),
        temp_dir=temp_dir,
        work_dir=work_dir,
        cache_dir=cache_dir,
        metadata=env_info,
    )


def is_containerized() -> bool:
    """判断是否在容器化环境中运行"""
    return detect_container_type() != ContainerType.NONE


def get_singularity_image_path(image_name: str) -> Optional[Path]:
    """
    获取 Singularity 镜像路径

    Args:
        image_name: 镜像名称或路径

    Returns:
        镜像路径或 None
    """
    if detect_container_type() != ContainerType.SINGULARITY:
        return None

    container_dir = os.environ.get("SINGULARITY_CONTAINER_DIR", "")

    if not container_dir:
        return None

    singularity_dir = Path(container_dir)

    possible_paths = [
        singularity_dir / image_name,
        singularity_dir / f"{image_name}.sif",
        singularity_dir / "images" / image_name,
        singularity_dir / "images" / f"{image_name}.sif",
    ]

    for path in possible_paths:
        if path.exists():
            return path

    return None
