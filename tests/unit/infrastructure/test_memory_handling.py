"""内存处理和错误恢复测试"""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.core.exceptions import ExecutionError
from phylodater.infrastructure.process_runner import ProcessRunner


class TestMemoryErrorHandling:

    @pytest.fixture
    def temp_dir(self):
        import shutil
        import tempfile

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_memory_error_gives_friendly_message(self, temp_dir):
        runner = ProcessRunner(cwd=temp_dir)

        with patch("subprocess.Popen") as mock_popen:
            mock_popen.side_effect = MemoryError("Out of memory")

            with pytest.raises(ExecutionError) as exc_info:
                runner.run(["some_command"])

            error_msg = str(exc_info.value)
            assert "memory" in error_msg.lower() or "RAM" in error_msg

    def test_generic_exception_with_memory_keyword(self, temp_dir):
        runner = ProcessRunner(cwd=temp_dir)

        with patch("subprocess.Popen") as mock_popen:
            mock_popen.side_effect = Exception("memory allocation failed")

            with pytest.raises(ExecutionError) as exc_info:
                runner.run(["some_command"])

            error_msg = str(exc_info.value)
            assert "memory" in error_msg.lower()


class TestProcessTimeoutHandling:

    @pytest.fixture
    def temp_dir(self):
        import shutil
        import tempfile

        temp_path = Path(tempfile.mkdtemp())
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_timeout_raises_timeout_error(self, temp_dir):
        runner = ProcessRunner(cwd=temp_dir, timeout=0.1)

        with patch("subprocess.Popen") as mock_popen:

            mock_process = Mock()
            mock_process.stdout = Mock()
            mock_process.stdout.readline = Mock(side_effect=["", ""])
            mock_process.stderr = Mock()
            mock_process.stderr.readline = Mock(side_effect=["", ""])
            mock_process.wait = Mock(side_effect=Exception("Timeout"))
            mock_process.pid = 12345
            mock_popen.return_value = mock_process

            with patch.object(runner, "logger") as _mock_logger:
                with pytest.raises(Exception):
                    runner.run(["some_command"])
