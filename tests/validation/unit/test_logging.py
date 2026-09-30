"""Unit tests for logging infrastructure."""

from phylodater.infrastructure import LogLevel, get_logger


class TestLogger:
    def test_singleton(self):
        logger1 = get_logger()
        logger2 = get_logger()
        assert logger1 is logger2

    def test_log_level_enum(self):
        assert LogLevel.DEBUG.value == "DEBUG"
        assert LogLevel.INFO.value == "INFO"

    def test_logger_level_set(self):
        logger = get_logger()
        logger.level = LogLevel.DEBUG
        assert logger.level == LogLevel.DEBUG

    def test_set_base_dir(self, tmp_test_dir):
        logger = get_logger()
        logger.set_base_dir(tmp_test_dir)
        assert logger._base_dir == tmp_test_dir
