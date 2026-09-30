"""ProcessRunner 输出体积护栏与返回码检查测试

覆盖审阅项：
- C-34：stdout/stderr 两个读线程共享的输出体积计数器必须线程安全
- C-35：``run(check=True)`` 在返回码不在 ``expected_returncodes`` 时抛错
"""

import inspect
import sys
import threading

import pytest

from phylodater.core.exceptions import ExecutionError
from phylodater.infrastructure.process_runner import ProcessRunner, _OutputBudget

# 每行 100 个字符；ProcessRunner 按 len(line) + 1 记账（换行符计入）
_LINE = "x" * 100
_LINES_PER_STREAM = 30
_BYTES_PER_STREAM = _LINES_PER_STREAM * (len(_LINE) + 1)  # 3030


def _emit_script(both_streams: bool) -> str:
    """生成一段向 stdout（可选同时向 stderr）刷定长输出的 python -c 脚本。"""
    script = (
        "import sys\n"
        f"line = {_LINE!r}\n"
        f"for _ in range({_LINES_PER_STREAM}):\n"
        "    print(line)\n"
    )
    if both_streams:
        script += (
            f"for _ in range({_LINES_PER_STREAM}):\n"
            "    print(line, file=sys.stderr)\n"
        )
    return script


@pytest.fixture
def temp_dir(tmp_path):
    return tmp_path


class TestOutputBudgetAccounting:
    """C-34：计数器必须不丢更新，且超限判定只发生一次。"""

    def test_no_updates_lost_under_concurrent_records(self):
        """8 线程各累加 2000 次，合并总量必须精确等于 16000。

        诚实标注：在 GIL 打开的 CPython 上这条断言**并不能证伪旧实现**——
        解释器把 LOAD/BINARY_OP/STORE 融合成 macro op，实测 8x3000 的裸 ``+=``
        也不丢更新（见 _OutputBudget 的文档）。本测试锁的是护栏计数的**契约**
        （每次累加都可见、总量精确），而不是复现一个只在自由线程构建（PEP 703）
        或别的解释器上才会显形的竞态。真正证伪"按线程各自记账"的是
        :meth:`test_output_limit_requires_both_streams`。
        """
        budget = _OutputBudget(max_size=10**9)
        n_threads, per_thread = 8, 2000
        barrier = threading.Barrier(n_threads)

        def worker():
            barrier.wait()
            for _ in range(per_thread):
                budget.record(1)

        # 尽量提高交错概率（即使调度不配合，断言依然成立）
        previous_interval = sys.getswitchinterval()
        sys.setswitchinterval(1e-6)
        try:
            threads = [threading.Thread(target=worker) for _ in range(n_threads)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
        finally:
            sys.setswitchinterval(previous_interval)

        assert budget.total == n_threads * per_thread
        assert not budget.exceeded

    def test_exceed_detected_exactly_once(self):
        """越线只被判定一次，且由累加后真正越过上限的那次 record 返回 True。"""
        limit = 100
        budget = _OutputBudget(max_size=limit)
        triggered = []
        lock = threading.Lock()

        def worker():
            for _ in range(20):
                if budget.record(1):
                    with lock:
                        triggered.append(True)

        threads = [threading.Thread(target=worker) for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # 8 * 20 = 160 次 +1，上限 100 → 必然越线，但只报一次
        assert budget.exceeded
        assert len(triggered) == 1
        assert budget.total == 160
        assert isinstance(budget.error, ExecutionError)
        assert str(limit) in str(budget.error)

    def test_output_limit_requires_both_streams(self, temp_dir):
        """护栏按两线程之和判定：单流不越线、双流越线。

        这是 C-34 的行为规格。若把上限按线程各自折算（或反过来完全丢掉另一
        线程的增量），第一或第二个断言就会失败。
        """
        single = ProcessRunner(
            cwd=temp_dir,
            monitor_memory=False,
            max_output_size=_BYTES_PER_STREAM * 2,  # 6060 > 单流 3030
        )
        result = single.run([sys.executable, "-c", _emit_script(both_streams=False)])
        # 只有一路输出时不应触发护栏（也证明上限没被"对半砍"）
        assert result.returncode == 0
        assert len(result.stdout.splitlines()) == _LINES_PER_STREAM

        both = ProcessRunner(
            cwd=temp_dir,
            monitor_memory=False,
            max_output_size=_BYTES_PER_STREAM + _BYTES_PER_STREAM // 2,  # 4545
        )
        with pytest.raises(ExecutionError) as exc_info:
            both.run([sys.executable, "-c", _emit_script(both_streams=True)])
        assert "Output size exceeded maximum limit" in str(exc_info.value)

    def test_generous_limit_keeps_all_output(self, temp_dir):
        """未越线时两路输出都完整保留（回归保护：加锁改动不得吞行）。"""
        runner = ProcessRunner(
            cwd=temp_dir, monitor_memory=False, max_output_size=1024 * 1024
        )
        result = runner.run([sys.executable, "-c", _emit_script(both_streams=True)])
        assert result.returncode == 0
        assert len(result.stdout.splitlines()) == _LINES_PER_STREAM
        assert len(result.stderr.splitlines()) == _LINES_PER_STREAM


class TestRunCheckReturnCode:
    """C-35：``check`` 参数把"必须记得检查 returncode"变成"默认安全"。"""

    def test_check_defaults_to_false_for_backward_compatibility(self):
        """默认 False 是刻意的：treePL prime/CV、MCMCTree 速率预估等调用点按设计
        需要拿到非零码与其输出；改成默认 True 会把它们的"降级路径"变成异常。
        """
        param = inspect.signature(ProcessRunner.run).parameters["check"]
        assert param.default is False

    def test_unexpected_returncode_returned_when_check_false(self, temp_dir):
        runner = ProcessRunner(cwd=temp_dir, monitor_memory=False)
        result = runner.run([sys.executable, "-c", "import sys; sys.exit(4)"])
        assert result.returncode == 4

    def test_unexpected_returncode_raises_when_check_true(self, temp_dir):
        runner = ProcessRunner(cwd=temp_dir, monitor_memory=False)
        script = (
            "import sys; "
            "print('boom: cannot open seed.dat', file=sys.stderr); "
            "sys.exit(4)"
        )
        cmd = [sys.executable, "-c", script]
        with pytest.raises(ExecutionError) as exc_info:
            runner.run(cmd, check=True)

        message = str(exc_info.value)
        assert "returncode 4" in message
        assert "boom: cannot open seed.dat" in message
        # 分类结论要随异常一起传出去，调用方不读日志也能拿到诊断依据
        assert "Suggestion:" in message

    def test_expected_returncode_tolerated_with_check_true(self, temp_dir):
        """PATHd8 式的"退出 1 仍算成功"与强制检查可以共存。"""
        runner = ProcessRunner(cwd=temp_dir, monitor_memory=False)
        result = runner.run(
            [sys.executable, "-c", "import sys; sys.exit(1)"],
            expected_returncodes={0, 1},
            check=True,
        )
        assert result.returncode == 1

    def test_zero_returncode_passes_with_check_true(self, temp_dir):
        runner = ProcessRunner(cwd=temp_dir, monitor_memory=False)
        result = runner.run([sys.executable, "-c", "print('ok')"], check=True)
        assert result.returncode == 0
        assert result.stdout.strip() == "ok"

    def test_check_raise_does_not_leak_the_child_process(self, temp_dir):
        """抛 ExecutionError 前已清理子进程（``_cleanup`` 会把句柄置空）。"""
        runner = ProcessRunner(cwd=temp_dir, monitor_memory=False)
        with pytest.raises(ExecutionError):
            runner.run([sys.executable, "-c", "import sys; sys.exit(7)"], check=True)
        assert runner._process is None
