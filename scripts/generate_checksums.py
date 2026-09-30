#!/usr/bin/env python3
"""
校验和生成工具

用于生成示例输出文件的校验和，供持续集成系统验证。

Usage:
    python scripts/generate_checksums.py <output_directory>
    python scripts/generate_checksums.py --recursive <directory>   # 写入 <directory>/checksums.sha256

Example:
    python scripts/generate_checksums.py examples/output
    python scripts/generate_checksums.py --recursive test-data/benchmark

``--recursive`` walks the tree, skips per-run output directories (``runs/``,
``__pycache__``) and the checksum file itself, and writes the result to
``<directory>/checksums.sha256`` so that ``sha256sum -c`` verifies from inside
that directory.
"""

import argparse
import hashlib
import sys
from pathlib import Path


def calculate_sha256(file_path: Path) -> str:
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def generate_checksums(output_dir: Path) -> dict:
    checksums = {}

    if not output_dir.exists():
        print(f"Error: Output directory not found: {output_dir}")
        return checksums

    for file_path in output_dir.iterdir():
        if file_path.is_file():
            checksum = calculate_sha256(file_path)
            checksums[file_path.name] = checksum
            print(f"{checksum}  {file_path.name}")

    return checksums


def update_checksum_file(checksums: dict, output_file: Path):
    lines = [
        "# PhyloDater Example Output Checksums",
        "# Generated for continuous integration validation",
        "#",
        "# Usage: sha256sum -c checksums.sha256",
        "#",
        f"# Generated: {__import__('datetime').datetime.now().isoformat()}",
        "",
        "# Expected output files from example run",
    ]

    for filename, checksum in sorted(checksums.items()):
        lines.append(f"{checksum}  {filename}")

    output_file.write_text("\n".join(lines))

    print(f"\nChecksums saved to: {output_file}")


def generate_recursive_checksums(root: Path) -> dict:
    """递归收集 ``root`` 下的文件校验和（键为相对 ``root`` 的路径）。

    跳过 ``runs/``（每次运行的产物，按需重新生成）、``__pycache__`` 与校验文件本身，
    否则清单会随一次普通运行不停变化，失去"验证数据集本身"的意义。
    """
    checksums = {}
    for file_path in sorted(root.rglob("*")):
        if not file_path.is_file():
            continue
        rel = file_path.relative_to(root).as_posix()
        parts = file_path.relative_to(root).parts
        if rel == "checksums.sha256" or "runs" in parts or "__pycache__" in parts:
            continue
        checksums[rel] = calculate_sha256(file_path)
    return checksums


def update_recursive_checksum_file(checksums: dict, root: Path) -> Path:
    lines = [
        "# PhyloDater dataset checksums",
        f"# Regenerate: python scripts/generate_checksums.py --recursive {root}",
        f"# Verify    : (cd {root} && sha256sum -c checksums.sha256)",
        "#",
        "# Per-run output under results/runs/ is deliberately excluded.",
    ]
    lines += [f"{digest}  {name}" for name, digest in sorted(checksums.items())]
    target = root / "checksums.sha256"
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("directory", type=Path, help="directory to checksum")
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="walk the tree and write <directory>/checksums.sha256",
    )
    args = parser.parse_args()

    root: Path = args.directory
    if not root.exists():
        print(f"Error: Output directory not found: {root}")
        sys.exit(1)

    print(f"Generating checksums for: {root}")
    print("-" * 70)

    if args.recursive:
        checksums = generate_recursive_checksums(root)
        for name, digest in sorted(checksums.items()):
            print(f"{digest}  {name}")
        target = update_recursive_checksum_file(checksums, root)
        print(f"\nChecksums saved to: {target}")
        print(f"\nTo verify:\n  (cd {root} && sha256sum -c checksums.sha256)")
        return

    checksums = generate_checksums(root)

    if not checksums:
        print("No files found in output directory.")
        sys.exit(1)

    checksum_file = Path("examples/checksums.sha256")
    update_checksum_file(checksums, checksum_file)

    print("\nTo verify outputs in CI:")
    print(f"  sha256sum -c {checksum_file}")


if __name__ == "__main__":
    main()
