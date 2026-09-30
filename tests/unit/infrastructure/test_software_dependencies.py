"""SoftwareDependencyManager 测试

覆盖审阅项：
- C-56：MD-Cat 的全名解释必须是"速率类别"，不是 Catmull-Friedman 样条
- C-57：依赖探测必须可移植（``shutil.which``），并区分"未安装"与"探测失败"
"""

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from phylodater.infrastructure import software_dependencies as sd
from phylodater.infrastructure.software_dependencies import (
    SOFTWARE_CATALOG,
    DependencyProbe,
    ProbeStatus,
    SoftwareDependencyManager,
)

# 一个不会是任何真实程序的名字；Windows 上旧实现连 `which` 都找不到，
# 于是所有依赖（包括这个不存在的）都同样报 False，用户无从分辨原因。
MISSING_NAME = "phylodater_definitely_not_an_installed_executable_zzz"
# sys.executable 是绝对路径，shutil.which 对含目录分隔符的入参直接做
# isfile + X_OK 检查，因此在任何平台都能被稳定"探测为已安装"。
KNOWN_EXECUTABLE = sys.executable


def _always_raise(exc):
    """供 monkeypatch 使用的"总是抛异常"替身。"""

    def _inner(*args, **kwargs):
        raise exc

    return _inner


class _RecordingLogger:
    """记录 warning 文本，用来验证"探测失败"被单独说出来了。"""

    def __init__(self):
        self.warnings = []

    def warning(self, message, *args, **kwargs):
        self.warnings.append(str(message))


class TestMDCatCatalogEntry:
    """C-56：目录解释会出现在自检/安装提示里，必须与上游语义一致。"""

    def test_description_is_rate_categories_not_spline(self):
        description = SOFTWARE_CATALOG["mdcat"].description
        assert "rate Categories" in description
        assert "categorical" in description.lower()
        # 旧描述把 MD-Cat 说成样条定年法（Catmull-Rom / Friedman 密集比样条）
        assert "Catmull" not in description
        assert "Friedman" not in description
        assert "spline" not in description.lower()

    def test_no_catalog_entry_reintroduces_the_wrong_spline_story(self):
        for info in SOFTWARE_CATALOG.values():
            assert "Catmull" not in info.description

    def test_upstream_readme_wording_is_the_basis_for_the_fix(self):
        """判据来自项目自带的上游 README（只读核对，文件不在则跳过）。

        上游目录 ``PhyloDater-参考软件`` 与被审目录 ``PhyloDater-项目代码`` 同级，
        故从本文件上溯 4 层到 ``PhyloDater/``。
        """
        reference = (
            Path(__file__).resolve().parents[4]
            / "PhyloDater-参考软件"
            / "原始代码-MD-Cat"
            / "README.md"
        )
        if not reference.exists():
            pytest.skip(f"上游 MD-Cat README 不在本检出中: {reference}")
        text = reference.read_text(encoding="utf-8")
        assert "categorical distribution" in text
        assert "number of rate categories" in text

    def test_executable_and_conda_name_mismatch_is_recorded_not_hidden(self):
        """记录（而非静默改）命名不一致：上游入口是 ``md_cat.py``。

        本目录的 ``executable``/``conda_package`` 都写作 ``mdcat``，而上游
        setup.py 以 ``scripts=['md_cat.py', 'simulate.py']`` 安装、发行名是
        ``MD-Cat``（见 MD_Cat.egg-info/PKG-INFO）；anaconda.org 上能查到的
        ``mdcat`` 是 conda-forge 的终端 Markdown 查看器，与 MD-Cat 无关。
        适配器实际探测/调用的路径来自配置项 ``software_paths.mdcat_bin``
        （代码默认 ``mdcat``，config.example.yaml 写 ``MD-Cat``）。三处对齐
        涉及 adapters/、infrastructure/configuration.py 与文档，属跨文件改动，
        见交付说明。本测试只钉住"当前值 + 已知不一致"这一事实，避免将来
        有人只改一处造成第四种写法。
        """
        info = SOFTWARE_CATALOG["mdcat"]
        assert info.executable == "mdcat"  # 已知：与上游入口 md_cat.py 不一致
        assert info.conda_package == "mdcat"
        assert "mdcat" in SoftwareDependencyManager.get_install_command("mdcat")

    def test_get_install_command_still_resolves_mdcAT(self):
        """改描述不得破坏按名字查目录（大小写/连字符/下划线归一化）。"""
        for query in ("mdcat", "MD-Cat", "MD_Cat", "md_cat"):
            command = SoftwareDependencyManager.get_install_command(query)
            assert "not found in catalog" not in command
            assert "mdcat" in command


class TestPortableDependencyProbe:
    """C-57：探测必须跨平台，且不把"探测失败"混同为"未安装"。"""

    def test_probe_reports_available_for_real_executable(self):
        probes = SoftwareDependencyManager.probe_dependencies([KNOWN_EXECUTABLE])
        probe = probes[KNOWN_EXECUTABLE]
        assert probe.status is ProbeStatus.AVAILABLE
        assert probe.available is True
        assert probe.path and Path(probe.path).exists()

    def test_probe_reports_missing_for_unknown_name(self):
        probes = SoftwareDependencyManager.probe_dependencies([MISSING_NAME])
        probe = probes[MISSING_NAME]
        assert probe.status is ProbeStatus.MISSING
        assert probe.available is False
        assert probe.error is None

    def test_probe_failure_is_not_reported_as_missing(self, monkeypatch):
        """PATH 不可读之类是"探测没做成"（可用性未知），不是"确认未装"。"""
        monkeypatch.setattr(
            sd.shutil,
            "which",
            _always_raise(PermissionError("Permission denied: '/nope'")),
        )
        probes = SoftwareDependencyManager.probe_dependencies([MISSING_NAME])
        probe = probes[MISSING_NAME]
        assert probe.status is ProbeStatus.PROBE_FAILED
        assert probe.available is False  # 未知也不说成"可用"（fail-closed 的一侧）
        assert "Permission denied" in probe.error
        assert "probe failed" in str(probe)

    def test_bool_view_collapsesto_false_but_says_why(self, monkeypatch):
        """check_all_dependencies 仍是 {软件: bool}，探测失败额外给一条 WARNING。"""
        recorder = _RecordingLogger()
        logging_module = importlib.import_module("phylodater.infrastructure.logging")
        monkeypatch.setattr(logging_module, "get_logger", lambda *a, **k: recorder)
        monkeypatch.setattr(
            sd.SoftwareDependencyManager,
            "probe_dependencies",
            staticmethod(
                lambda names: {
                    n: DependencyProbe(n, ProbeStatus.PROBE_FAILED, error="boom")
                    for n in names
                }
            ),
        )

        results = SoftwareDependencyManager.check_all_dependencies(["treePL"])

        assert results == {"treePL": False}
        assert len(recorder.warnings) == 1
        assert "probe failed" in recorder.warnings[0].lower()
        assert "UNKNOWN" in recorder.warnings[0]
        assert "treePL" in recorder.warnings[0]

    def test_probe_does_not_shell_out_to_which(self, monkeypatch):
        """回归护栏：不得再用 POSIX 外部命令 ``which`` 探测（Windows 上必炸）。"""
        refuse = _always_raise(
            AssertionError("依赖探测不应起子进程（shutil.which 才是可移植正解）")
        )
        monkeypatch.setattr(subprocess, "run", refuse)
        monkeypatch.setattr(subprocess, "Popen", refuse)

        results = SoftwareDependencyManager.check_all_dependencies(
            [MISSING_NAME, KNOWN_EXECUTABLE]
        )
        assert results[MISSING_NAME] is False
        assert results[KNOWN_EXECUTABLE] is True

    def test_bool_view_matches_shutil_which(self):
        names = [MISSING_NAME, KNOWN_EXECUTABLE]
        results = SoftwareDependencyManager.check_all_dependencies(names)
        for name in names:
            assert results[name] == (shutil.which(name) is not None)

    def test_empty_input_returns_empty_mapping(self):
        assert SoftwareDependencyManager.check_all_dependencies([]) == {}
        assert SoftwareDependencyManager.probe_dependencies([]) == {}

    def test_probe_helpers_are_reachable_without_subprocess(self):
        """探测走的是 shutil.which：模块级的 subprocess 依赖应当不存在。"""
        assert sd.shutil.which is shutil.which
        assert "subprocess" not in getattr(sd, "__dict__", {})
