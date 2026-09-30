"""Functional tests for the 'phylodater check' CLI."""

from tests.validation.helpers import run_cli


class TestCliCheck:
    def test_check_runs(self):
        result = run_cli(["check"])
        assert result.returncode == 0, result.stderr

    def test_self_test_alias(self):
        result = run_cli(["self-test"])
        assert result.returncode == 0, result.stderr
