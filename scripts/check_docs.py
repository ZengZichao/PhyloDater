#!/usr/bin/env python3
"""Check the repository's documentation policy (links, bilingual pairs, banners).

Enforces what docs/README.md promises:

* every internal markdown link resolves to a file that exists;
* every document has both an English (bare name, canonical) and a Chinese
  (``_CN``) version, and no ``*_EN.md`` duplicate comes back;
* the language switch line at the top of a document lists English first.

Documented exemptions (see the policy table in docs/README.md):

* ``.github/ISSUE_TEMPLATE/*.md`` — one file per picker entry, bilingual inline.

Exit code 0 = policy satisfied, 1 = violations printed.

Usage:  python scripts/check_docs.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LINK_RE = re.compile(r"\[[^\]]*\]\(([^)\s]+)\)")
EXEMPT_PAIRING = (".github/ISSUE_TEMPLATE/", ".pytest_cache/", "build/")


def markdown_files() -> list[Path]:
    return sorted(
        p
        for p in REPO_ROOT.rglob("*.md")
        if ".git/" not in str(p.relative_to(REPO_ROOT)) and "__pycache__" not in str(p)
    )


def check_links(files: list[Path]) -> list[str]:
    problems = []
    for md in files:
        for target in LINK_RE.findall(md.read_text(encoding="utf-8")):
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            relative = target.split("#")[0].strip()
            if not relative:
                continue
            if not (md.parent / relative).resolve().exists():
                problems.append(f"{md.relative_to(REPO_ROOT)} -> {target}")
    return problems


def check_pairing(files: list[Path]) -> list[str]:
    problems = []
    for md in files:
        relative = md.relative_to(REPO_ROOT).as_posix()
        if any(relative.startswith(prefix) for prefix in EXEMPT_PAIRING):
            continue
        stem = md.name[:-3]
        if stem.endswith("_EN"):
            problems.append(
                f"{relative}: duplicate English variant, fold it into '{stem[:-3]}.md'"
            )
            continue
        twin = md.parent / (
            (stem[:-3] if stem.endswith("_CN") else stem) + "_CN.md"
            if not stem.endswith("_CN")
            else stem[:-3] + ".md"
        )
        if not twin.exists():
            language = "English" if stem.endswith("_CN") else "Chinese"
            problems.append(f"{relative}: missing {language} twin ({twin.name})")
    return problems


def check_banners(files: list[Path]) -> list[str]:
    problems = []
    for md in files:
        relative = md.relative_to(REPO_ROOT).as_posix()
        lines = [
            line
            for line in md.read_text(encoding="utf-8").lstrip().splitlines()
            if line.strip() and not line.startswith("<!--")
        ]
        if not lines:
            continue
        head = lines[0]
        if "Language" not in head and "语言" not in head:
            continue
        english, chinese = head.find("English"), head.find("中文")
        if chinese != -1 and (english == -1 or chinese < english):
            problems.append(f"{relative}: language switch does not list English first")
    return problems


def main() -> int:
    files = markdown_files()
    sections = {
        "broken internal links": check_links(files),
        "bilingual pairing problems": check_pairing(files),
        "language switches not English-first": check_banners(files),
    }
    total = 0
    for title, problems in sections.items():
        print(f"{title}: {len(problems)}")
        for problem in problems:
            print(f"  {problem}")
        total += len(problems)
    if total:
        print(f"\ndocumentation policy violated in {total} place(s)")
        return 1
    print(f"\ndocumentation policy satisfied ({len(files)} markdown files checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
