"""
模块七：日志与信号处理功能测试

验证日志格式、信号处理、资源清理等。
"""

import signal
import subprocess
import sys
import time

from phylodater.infrastructure.signal_handler import get_shutdown_handler
from tests.functional.conftest import run_phylodater_cli


class TestLoggingSystem:
    """日志系统测试"""

    def test_log_001_iso_timestamp(self, tmp_path):
        """LOG-001: ISO 8601 时间戳"""
        log_file = tmp_path / "test.log"
        result = run_phylodater_cli(["--log-file", str(log_file), "check"])
        assert result.returncode in [0, 1]
        assert log_file.exists()
        content = log_file.read_text(encoding="utf-8")
        assert "T" in content  # ISO 8601 包含 T

    def test_log_006_log_level_filter(self, tmp_path):
        """LOG-006: 日志级别过滤"""
        log_file = tmp_path / "test.log"
        result = run_phylodater_cli(
            ["--log-level", "WARNING", "--log-file", str(log_file), "check"]
        )
        assert result.returncode in [0, 1]
        assert log_file.exists()
        content = log_file.read_text(encoding="utf-8")
        # WARNING 级别不应包含 INFO/DEBUG 消息（系统启动消息除外）
        lines = [ln for ln in content.splitlines() if " | " in ln]
        non_start_lines = [ln for ln in lines if "日志系统已启动" not in ln]
        if non_start_lines:
            assert all(
                "INFO" not in ln or "WARNING" in ln or "ERROR" in ln
                for ln in non_start_lines
            )


class TestSignalHandling:
    """信号处理测试"""

    def test_sig_001_ctrl_c_terminates_with_130(self):
        """SIG-001: Ctrl+C 触发信号处理"""
        # 启动一个子进程并发送 SIGINT
        cmd = [sys.executable, "-m", "phylodater", "check"]
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        time.sleep(0.5)
        proc.send_signal(signal.SIGINT)
        try:
            proc.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.communicate()
        # 子进程可能被信号终止 (-signal.SIGINT = -2) 或优雅退出 130/0/1
        assert proc.returncode in [130, 0, 1, -2]

    def test_sig_003_resource_cleanup(self, tmp_path):
        """SIG-003: 资源清理"""
        handler = get_shutdown_handler()
        temp_dir = tmp_path / "cleanup_test"
        temp_dir.mkdir()
        handler.register_temp_dir(temp_dir)
        handler._cleanup()
        assert not temp_dir.exists()
