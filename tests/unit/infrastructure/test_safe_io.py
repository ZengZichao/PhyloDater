"""``infrastructure.safe_io`` 的单元测试：路径校验与打开语义。"""

from __future__ import annotations

from pathlib import Path

import pytest

from phylodater.infrastructure.safe_io import safe_writer


def test_write_creates_file_with_expected_content(tmp_path: Path) -> None:
    target = tmp_path / "out.txt"
    with safe_writer(target, encoding="utf-8") as f:
        f.write("内容\n")
    assert target.read_text(encoding="utf-8") == "内容\n"


def test_rejects_parent_segment(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"\.\."):
        safe_writer(tmp_path / ".." / "escape.txt")


def test_rejects_relative_parent_segment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match=r"\.\."):
        safe_writer(Path("a") / ".." / "escape.txt")


def test_binary_mode(tmp_path: Path) -> None:
    target = tmp_path / "out.bin"
    with safe_writer(target, "wb") as f:
        f.write(b"\x00\x01")
    assert target.read_bytes() == b"\x00\x01"


def test_kwargs_passthrough_newline(tmp_path: Path) -> None:
    target = tmp_path / "out.tsv"
    with safe_writer(target, newline="", encoding="utf-8") as f:
        f.write("a\r\nb\n")
    assert b"\r\n" in target.read_bytes()


def test_default_mode_is_text(tmp_path: Path) -> None:
    target = tmp_path / "plain.txt"
    with safe_writer(target) as f:
        f.write("x")
    assert target.read_text(encoding="utf-8") == "x"
