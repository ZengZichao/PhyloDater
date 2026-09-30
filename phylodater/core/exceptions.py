"""
PhyloDater 异常层次结构

所有异常继承自 PhyloDaterError，形成清晰的错误分类体系。
错误分类支持可操作的错误信息，便于用户快速定位和解决问题。
"""

from enum import Enum
from typing import Optional


class ErrorCategory(Enum):
    """错误分类枚举

    用于区分不同性质的错误，帮助用户和系统针对性处理。
    """

    USER_INPUT = "user_input"  # 用户输入错误（文件格式、参数错误等）
    ENGINE_FAILURE = "engine_failure"  # 引擎运行失败（外部软件崩溃、收敛失败等）
    RESOURCE_LIMITED = "resource_limited"  # 系统资源不足（内存、超时等）
    DATA_ISSUE = "data_issue"  # 数据问题（数据损坏、不一致等）
    SYSTEM_ERROR = "system_error"  # 系统级错误（文件权限、环境配置等）
    UNKNOWN = "unknown"  # 未知错误


class PhyloDaterError(Exception):
    """所有 PhyloDater 异常的基类"""

    error_category: ErrorCategory = ErrorCategory.UNKNOWN
    recoverable: bool = False
    suggestion: Optional[str] = None

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        super().__init__(message)
        self.suggestion = suggestion

    def to_dict(self) -> dict:
        """转换为字典，用于结构化日志或序列化"""
        return {
            "type": self.__class__.__name__,
            "category": self.error_category.value,
            "message": str(self),
            "recoverable": self.recoverable,
            "suggestion": self.suggestion,
        }


# ==================== 验证错误 ====================


class ValidationError(PhyloDaterError):
    """输入验证失败"""

    error_category = ErrorCategory.USER_INPUT
    recoverable = False

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请检查输入文件格式是否正确，参考文档中的输入要求。"
        super().__init__(message, suggestion)


class TreeValidationError(ValidationError):
    """树结构问题（无根树、缺失分支长度等）"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请确保树文件为有效的 Newick 格式，包含分支长度，且已定根。"
        super().__init__(message, suggestion)


class AlignmentValidationError(ValidationError):
    """比对问题（长度不一致、格式错误等）"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请确保比对文件为有效的 FASTA 格式，所有序列长度一致。"
        super().__init__(message, suggestion)


class ConsistencyError(ValidationError):
    """树与比对不一致"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请确保树文件中的物种名称与比对文件中的序列名称一致。"
        super().__init__(message, suggestion)


# ==================== 校准错误 ====================


class CalibrationError(PhyloDaterError):
    """校准点相关错误"""

    error_category = ErrorCategory.USER_INPUT


class ConfigurationError(PhyloDaterError):
    """配置错误"""

    error_category = ErrorCategory.USER_INPUT


class CalibrationConflictError(CalibrationError):
    """祖先-后代时序冲突"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请检查校准点的时间约束，确保祖先节点的最小年龄不超过后代节点的最大年龄。"
        super().__init__(message, suggestion)


class CalibrationResolutionError(CalibrationError):
    """MRCA 定位失败"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请确保校准点中指定的分类单元名称与树中的叶节点名称匹配。"
        super().__init__(message, suggestion)


# ==================== 执行错误 ====================


class ExecutionError(PhyloDaterError):
    """外部程序执行失败"""

    error_category = ErrorCategory.ENGINE_FAILURE
    recoverable = False

    def __init__(
        self,
        message: str = "",
        category: ErrorCategory = ErrorCategory.ENGINE_FAILURE,
        recoverable: bool = False,
        suggestion: Optional[str] = None,
    ) -> None:
        self.error_category = category
        self.recoverable = recoverable
        if suggestion is None:
            suggestion = self._default_suggestion(category)
        super().__init__(message, suggestion)

    @staticmethod
    def _default_suggestion(category: ErrorCategory) -> str:
        """根据错误类别返回默认建议"""
        suggestions = {
            ErrorCategory.USER_INPUT: "请检查输入参数是否正确。",
            ErrorCategory.ENGINE_FAILURE: "请检查外部软件是否正确安装，或尝试使用其他定年方法。",
            ErrorCategory.RESOURCE_LIMITED: "请尝试减少数据规模或增加系统资源（内存/CPU）。",
            ErrorCategory.DATA_ISSUE: "请检查数据完整性，必要时重新下载或生成数据。",
            ErrorCategory.SYSTEM_ERROR: "请检查系统环境配置和文件权限。",
            ErrorCategory.UNKNOWN: "请查看详细日志或联系开发者报告此问题。",
        }
        return suggestions.get(category, "请查看详细日志获取更多信息。")


class ExecutionTimeoutError(ExecutionError):
    """超时终止"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = (
                "请尝试增加超时时间（环境变量 PHYLODATER_TIMEOUT），或减少数据规模。"
            )
        super().__init__(
            message,
            category=ErrorCategory.RESOURCE_LIMITED,
            recoverable=True,
            suggestion=suggestion,
        )


class ConvergenceWarning(ExecutionError):
    """MCMC 未充分收敛（非致命）"""

    error_category = ErrorCategory.ENGINE_FAILURE
    recoverable = True

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请尝试增加 MCMC 迭代次数（burnin/sampfreq/nsample 参数）或增加运行次数。"
        super().__init__(
            message,
            category=ErrorCategory.ENGINE_FAILURE,
            recoverable=True,
            suggestion=suggestion,
        )


# ==================== 结果错误 ====================


class ResultParsingError(PhyloDaterError):
    """结果解析失败"""

    error_category = ErrorCategory.ENGINE_FAILURE
    recoverable = False

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "外部软件可能未正常完成，请检查日志文件和原始输出。"
        super().__init__(message, suggestion)


# ==================== 警告类（非致命） ====================


class SemanticDegradationWarning(PhyloDaterError, UserWarning):
    """约束语义降级（非致命，收集到摘要）

    双重继承：既是 Exception（可通过 raise 捕获）又是 UserWarning（可通过 warnings.warn 发出）。
    适配器中使用 warnings.warn() 发出此警告，约束降级值仍然正常返回。
    """

    error_category = ErrorCategory.USER_INPUT
    recoverable = True

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = (
                "部分校准约束被降级处理，结果可能存在偏差。建议检查校准点配置。"
            )
        super().__init__(message, suggestion)


class EnvironmentWarning(PhyloDaterError):
    """软件版本未验证或不兼容（非致命）"""

    error_category = ErrorCategory.SYSTEM_ERROR
    recoverable = True

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请确保所有外部软件已正确安装并添加到系统路径。"
        super().__init__(message, suggestion)


class ResourceWarning(PhyloDaterError):
    """资源使用超阈值（非致命）"""

    error_category = ErrorCategory.RESOURCE_LIMITED
    recoverable = True

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请监控资源使用情况，必要时减少数据规模或增加系统资源。"
        super().__init__(message, suggestion)


class UnknownParameterWarning(PhyloDaterError):
    """未知参数警告"""

    error_category = ErrorCategory.USER_INPUT
    recoverable = True

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请检查参数拼写是否正确，或查看文档获取支持的参数列表。"
        super().__init__(message, suggestion)


# ==================== 其他错误 ====================


class UnknownMethodError(PhyloDaterError):
    """DatingMethodRegistry 中未找到指定方法名"""

    error_category = ErrorCategory.USER_INPUT
    recoverable = False

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请检查方法名称拼写是否正确。支持的定年方法请参考文档。"
        super().__init__(message, suggestion)


class PhyloFormatError(ValidationError):
    """不支持的文件格式或格式错误"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请将文件转换为支持的格式。FASTA格式为推荐格式。"
        super().__init__(message, suggestion)


class TaxonomyConflictError(CalibrationError):
    """分类学信息冲突（标签不一致、循环依赖等）"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "请检查分类学表中的名称与树中的末端标签是否一致。"
        super().__init__(message, suggestion)


class MonophylyError(CalibrationError):
    """单系群校验失败"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = (
                "指定的分类群在树中不构成单系群，请检查校准点定义或树拓扑结构。"
            )
        super().__init__(message, suggestion)


class NegativeBranchLengthError(TreeValidationError):
    """负分支长度错误（退出码3）"""

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "树文件包含负分支长度，请检查树文件或使用 --skip-length-check 跳过检查。"
        super().__init__(message, suggestion)


class InputSizeLimitError(ValidationError):
    """输入文件大小超限"""

    error_category = ErrorCategory.RESOURCE_LIMITED
    recoverable = True

    def __init__(self, message: str = "", suggestion: Optional[str] = None) -> None:
        if suggestion is None:
            suggestion = "输入文件过大，请减小数据规模或使用 --low-memory 模式。"
        super().__init__(message, suggestion)


# ==================== 错误分类辅助函数 ====================


def classify_execution_error(stderr: str, returncode: int) -> tuple:
    """根据外部程序输出版本分类执行错误

    Args:
        stderr: 标准错误输出
        returncode: 进程返回码

    Returns:
        (ErrorCategory, str): 错误类别和建议消息
    """
    stderr_lower = stderr.lower()

    if returncode == -9 or returncode == 137:
        return (
            ErrorCategory.RESOURCE_LIMITED,
            "进程被系统强制终止（可能是内存不足）。建议增加可用内存或减少数据规模。",
        )

    if "out of memory" in stderr_lower or "memoryerror" in stderr_lower:
        return (
            ErrorCategory.RESOURCE_LIMITED,
            "内存不足。请尝试减少数据规模、增加交换空间或使用更大内存的机器。",
        )

    if "timeout" in stderr_lower or "timed out" in stderr_lower:
        return (
            ErrorCategory.RESOURCE_LIMITED,
            "执行超时。请尝试增加超时时间或减少数据规模。",
        )

    if returncode == 127:
        return (
            ErrorCategory.SYSTEM_ERROR,
            "外部程序未找到。请确保软件已安装并添加到系统路径。",
        )

    if "permission denied" in stderr_lower:
        return (
            ErrorCategory.SYSTEM_ERROR,
            "权限不足。请检查文件权限或以更高权限运行。",
        )

    if "file not found" in stderr_lower or "no such file" in stderr_lower:
        return (
            ErrorCategory.USER_INPUT,
            "必要的输入文件不存在。请检查文件路径是否正确。",
        )

    if "error" in stderr_lower and "parse" in stderr_lower:
        return (
            ErrorCategory.USER_INPUT,
            "输入文件格式错误。请检查文件是否符合要求的格式。",
        )

    return (
        ErrorCategory.ENGINE_FAILURE,
        "外部程序执行失败。请查看详细日志或联系开发者。",
    )
