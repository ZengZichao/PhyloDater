"""
异常类的单元测试
"""

import pytest

from phylodater.core.exceptions import (
    CalibrationError,
    ConfigurationError,
    ExecutionError,
    PhyloDaterError,
    SemanticDegradationWarning,
    TreeValidationError,
    UnknownMethodError,
    ValidationError,
)


class TestPhyloDaterError:
    """测试基础异常类"""

    def test_basic_exception(self):
        """测试基本异常"""
        with pytest.raises(PhyloDaterError):
            raise PhyloDaterError("Test error")

    def test_exception_message(self):
        """测试异常消息"""
        error = PhyloDaterError("Custom message")
        assert str(error) == "Custom message"

    def test_exception_inheritance(self):
        """测试继承自 Exception"""
        assert issubclass(PhyloDaterError, Exception)


class TestValidationError:
    """测试验证错误"""

    def test_is_phylo_dater_error(self):
        """测试是 PhyloDaterError 的子类"""
        assert issubclass(ValidationError, PhyloDaterError)

    def test_raise_validation_error(self):
        """测试抛出验证错误"""
        with pytest.raises(ValidationError):
            raise ValidationError("Validation failed")


class TestTreeValidationError:
    """测试树验证错误"""

    def test_is_validation_error(self):
        """测试是 ValidationError 的子类"""
        assert issubclass(TreeValidationError, ValidationError)

    def test_raise_tree_error(self):
        """测试抛出树验证错误"""
        with pytest.raises(TreeValidationError):
            raise TreeValidationError("Invalid tree structure")


class TestCalibrationError:
    """测试校准错误"""

    def test_is_phylo_dater_error(self):
        """测试是 PhyloDaterError 的子类"""
        assert issubclass(CalibrationError, PhyloDaterError)

    def test_raise_calibration_error(self):
        """测试抛出校准错误"""
        with pytest.raises(CalibrationError):
            raise CalibrationError("Invalid calibration point")


class TestExecutionError:
    """测试执行错误"""

    def test_is_phylo_dater_error(self):
        """测试是 PhyloDaterError 的子类"""
        assert issubclass(ExecutionError, PhyloDaterError)

    def test_raise_execution_error(self):
        """测试抛出执行错误"""
        with pytest.raises(ExecutionError):
            raise ExecutionError("Command failed")


class TestUnknownMethodError:
    """测试未知方法错误"""

    def test_is_phylo_dater_error(self):
        """测试是 PhyloDaterError 的子类"""
        assert issubclass(UnknownMethodError, PhyloDaterError)

    def test_raise_unknown_method_error(self):
        """测试抛出未知方法错误"""
        with pytest.raises(UnknownMethodError):
            raise UnknownMethodError("Unknown method: test")


class TestConfigurationError:
    """测试配置错误"""

    def test_is_phylo_dater_error(self):
        """测试是 PhyloDaterError 的子类"""
        assert issubclass(ConfigurationError, PhyloDaterError)

    def test_raise_configuration_error(self):
        """测试抛出配置错误"""
        with pytest.raises(ConfigurationError):
            raise ConfigurationError("Invalid configuration")


class TestSemanticDegradationWarning:
    """测试语义降级警告"""

    def test_is_warning(self):
        """测试是 Warning 的子类"""
        assert issubclass(SemanticDegradationWarning, Warning)

    def test_raise_warning(self):
        """测试发出警告"""
        with pytest.warns(SemanticDegradationWarning):
            import warnings

            warnings.warn("Feature degraded", SemanticDegradationWarning)
