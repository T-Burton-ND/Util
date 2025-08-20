#!/usr/bin/env python3
"""
update_readme.py — Help-first README auto-updater
───────────────────────────────────────────────────────────────────────
Auto-refreshes README sections by scanning the repository for scripts and
capturing ONLY their `--help` / `-h` output to describe usage.

What it does
• Recursively discovers .py/.sh scripts, skipping junk (e.g., .git, venvs)
• For each script, runs:
    - Python:  python <script> -h
    - Bash:    <script> --help
• Extracts:
    – Title (filename), relative path
    – Language (Python/Bash)
    – Description = first non-empty line of help (before any "Usage" block)
    – Usage      = help text's "Usage" block (if present), otherwise empty
    – Lines of code, last modified time, category (by folder hint)
• Writes:
    – README.md blocks between markers:
        <!-- BEGIN AUTO-OVERVIEW --> ... <!-- END AUTO-OVERVIEW -->
        <!-- BEGIN AUTO-SCRIPTS  --> ... <!-- END AUTO-SCRIPTS  -->
    – Optional index files: script_index.json and script_index.csv
• Safe:
    – Dry-run by default (prints preview to stdout)
    – --write makes an in-place update with a .bak backup
    – Does not touch anything outside the markers

Usage
  python update_readme.py               # preview only
  python update_readme.py --write       # apply changes to README.md
  python update_readme.py --write-index # emit script_index.json/.csv

Recommended exclusions for speed/sanity:
  --exclude-dirs .git .venv build dist __pycache__ node_modules .mypy_cache
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Dict

# ─────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────
READ_ME = Path("README.md")
INCLUDE_EXT = {".py", ".sh"}
DEFAULT_EXCLUDE_DIRS = {
    ".git", ".svn", ".hg", ".mypy_cache", "__pycache__",
    ".venv", "venv", "env", "build", "dist", ".ipynb_checkpoints",
    ".pytest_cache", ".tox", ".idea", ".vscode", ".DS_Store",
}
IGNORE_HIDDEN = True
HELP_TIMEOUT = 7  # seconds per script

MARK_OVERVIEW_BEGIN = "<!-- BEGIN AUTO-OVERVIEW -->"
MARK_OVERVIEW_END   = "<!-- END AUTO-OVERVIEW -->"
MARK_SCRIPTS_BEGIN  = "<!-- BEGIN AUTO-SCRIPTS -->"
MARK_SCRIPTS_END    = "<!-- END AUTO-SCRIPTS -->"

CATEGORY_HINTS = {
    "automation": "Automation",
    "scraping": "Scraping",
    "visualization_scripts": "Visualization",
    "visualization": "Visualization",
    "viz": "Visualization",
}

# ─────────────────────────────────────────────────────────────────────
# Data model
# ─────────────────────────────────────────────────────────────────────
@dataclass
class ScriptInfo:
    path: Path
    language: str            # "Python" or "Bash"
    description: str
    usage: str
    loc: int
    modified: str            # ISO date string
    category: str

    @property
    def relpath(self) -> str:
        return str(self.path.as_posix())

    @property
    def name(self) -> str:
        return self.path.name

# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────
def is_hidden(p: Path) -> bool:
    return any(part.startswith(".") for part in p.parts)

def discover_scripts(root: Path, include_ext: set[str], exclude_dirs: set[str], ignore_hidden: bool = True) -> List[Path]:
    paths: List[Path] = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if ignore_hidden and is_hidden(p.relative_to(root)):
            continue
        if p.suffix.lower() in include_ext:
            if any(part in exclude_dirs for part in p.parts):
                continue
            paths.append(p)
    return sorted(paths)

def read_text(p: Path) -> str:
    try:
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""

def count_loc(text: str) -> int:
    return sum(1 for line in text.splitlines() if line.strip())

def last_modified_iso(p: Path) -> str:
    try:
        ts = p.stat().st_mtime
    except FileNotFoundError:
        return ""
    return dt.datetime.fromtimestamp(ts).strftime("%Y-%m-%d")

def extract_description_from_help(help_text: str) -> str:
    """
    Return the first non-empty line occurring before a 'Usage' header.
    """
    lines = [ln.strip("` ").rstrip() for ln in help_text.splitlines()]
    desc_lines = []
    for ln in lines:
        if not ln:
            # allow blank lines until we capture something meaningful
            if desc_lines:
                break
            continue
        if re.match(r"(?i)usage\s*:?", ln):
            break
        # skip obvious flag-only lines
        if ln.startswith("-"):
            continue
        desc_lines.append(ln)
        # keep it short: first non-empty human line is enough
        break
    return desc_lines[0] if desc_lines else "No --help detected"

def extract_usage_block(help_text: str) -> str:
    """
    Grab a compact 'Usage' block from the help text (first paragraph after 'Usage').
    """
    m = re.search(r"(?im)^\s*usage\s*:?\s*(.+)$", help_text)
    if not m:
        return ""
    # collect the 'Usage' line + any immediately following lines until a blank line
    start = m.start()
    tail = help_text[start:]
    block = []
    for ln in tail.splitlines():
        if not ln.strip():
            break
        block.append(ln.strip("` "))
        # avoid giant blocks
        if len(block) > 8:
            break
    # keep lines reasonably short in the table
    block = [ln if len(ln) <= 160 else (ln[:157] + "…") for ln in block]
    return "\n".join(block)

def try_help(path: Path) -> str:
    """
    Attempt to capture --help/-h. Python uses '-h', Bash uses '--help'.
    """
    try_cmds = []
    if path.suffix == ".py":
        try_cmds = [
            [sys.executable, str(path), "-h"],
            [sys.executable, str(path), "--help"],
        ]
    else:
        try_cmds = [
            [str(path), "--help"],
            [str(path), "-h"],
        ]

    for cmd in try_cmds:
        try:
            proc = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=HELP_TIMEOUT, check=False, text=True
            )
            out = (proc.stdout or "").strip()
            # accept non-empty output as 'help'
            if out:
                return out
        except Exception:
            continue
    return ""

def infer_category(path: Path) -> str:
    parts = [p.lower() for p in path.parts]
    for key, label in CATEGORY_HINTS.items():
        if key in parts:
            return label
    return "Other"

def scripts_to_markdown_overview(scripts: List[ScriptInfo]) -> str:
    n = len(scripts)
    n_py = sum(1 for s in scripts if s.language == "Python")
    n_sh = sum(1 for s in scripts if s.language == "Bash")
    today = dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    cats: Dict[str, int] = {}
    for s in scripts:
        cats[s.category] = cats.get(s.category, 0) + 1

    cats_line = " • ".join(f"{k}: {v}" for k, v in sorted(cats.items()))
    return (
        f"**Inventory summary** — {n} scripts "
        f"(Python: {n_py}, Bash: {n_sh}) • {cats_line} • Last scan: {today}\n"
    )

def scripts_to_markdown_tables(scripts: List[ScriptInfo]) -> str:
    # group by category
    by_cat: Dict[str, List[ScriptInfo]] = {}
    for s in scripts:
        by_cat.setdefault(s.category, []).append(s)
    out = []
    for cat in sorted(by_cat):
        out.append(f"### 🔹 {cat}\n")
        out.append("| Script | Lang | Description | Usage | LOC | Modified |")
        out.append("|---|---|---|---|---|---|")
        for s in sorted(by_cat[cat], key=lambda x: x.name.lower()):
            link = f"[`{s.name}`]({s.relpath})"
            desc = (s.description or "").replace("|", "\\|")
            usage = (s.usage or "").replace("|", "\\|")
            out.append(
                f"| {link} | {s.language} | {desc} | {usage} | {s.loc} | {s.modified} |"
            )
        out.append("")  # blank line between categories
    return "\n".join(out).strip()

def replace_block(text: str, begin: str, end: str, replacement: str) -> str:
    if begin not in text or end not in text:
        # If markers are missing, append a fresh block at the end.
        return text.rstrip() + f"\n\n{begin}\n{replacement}\n{end}\n"
    pattern = re.compile(
        re.escape(begin) + r"(.*?)" + re.escape(end),
        flags=re.DOTALL
    )
    return pattern.sub(begin + "\n" + replacement + "\n" + end, text)

def load_readme(path: Path) -> str:
    if path.exists():
        return read_text(path)
    # Create a minimal template if missing
    return f"# Util Repository – README\n\n## 📜 Script Overview\n\n{MARK_OVERVIEW_BEGIN}\n{MARK_OVERVIEW_END}\n\n## 🔧 Detailed Inventory\n\n{MARK_SCRIPTS_BEGIN}\n{MARK_SCRIPTS_END}\n"

def write_index_files(scripts: List[ScriptInfo], root: Path) -> None:
    if not scripts:
        return
    data = [asdict(s) for s in scripts]
    (root / "script_index.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    with (root / "script_index.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
        w.writeheader()
        for row in data:
            w.writerow(row)

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(description="Auto-update README tables from repository scripts (help-only).")
    ap.add_argument("--repo-root", type=Path, default=Path("."), help="Repository root (default: .)")
    ap.add_argument("--readme", type=Path, default=READ_ME, help="README path (default: README.md)")
    ap.add_argument("--include-ext", nargs="*", default=list(INCLUDE_EXT), help="File extensions to include")
    ap.add_argument("--exclude-dirs", nargs="*", default=list(DEFAULT_EXCLUDE_DIRS), help="Directory names to exclude")
    ap.add_argument("--no-ignore-hidden", action="store_true", help="Do not skip hidden files/folders")
    ap.add_argument("--write", action="store_true", help="Apply changes to README (otherwise print preview)")
    ap.add_argument("--write-index", action="store_true", help="Write script_index.json and script_index.csv")
    args = ap.parse_args()

    root: Path = args.repo_root.resolve()
    readme_path: Path = (root / args.readme).resolve()

    include_ext = {ext if ext.startswith(".") else f".{ext}" for ext in args.include_ext}
    exclude_dirs = set(args.exclude_dirs)
    ignore_hidden = not args.no_ignore_hidden

    # Discover scripts
    files = discover_scripts(root, include_ext, exclude_dirs, ignore_hidden)
    if not files:
        print("No scripts found. Adjust filters?")
        sys.exit(0)

    all_scripts: List[ScriptInfo] = []
    for path in files:
        text = read_text(path)                    # only for LOC
        lang = "Python" if path.suffix == ".py" else "Bash"

        help_text = try_help(path)
        description = extract_description_from_help(help_text) if help_text else "No --help detected"
        usage = extract_usage_block(help_text) if help_text else ""

        info = ScriptInfo(
            path=path.relative_to(root),
            language=lang,
            description=description,
            usage=usage,
            loc=count_loc(text),
            modified=last_modified_iso(path),
            category=infer_category(path),
        )
        all_scripts.append(info)

    # Build blocks
    overview_md = scripts_to_markdown_overview(all_scripts)
    tables_md   = scripts_to_markdown_tables(all_scripts)

    if args.write_index:
        write_index_files(all_scripts, root)

    # Load or make README; replace blocks
    content = load_readme(readme_path)
    new_content = replace_block(content, MARK_OVERVIEW_BEGIN, MARK_OVERVIEW_END, overview_md)
    new_content = replace_block(new_content, MARK_SCRIPTS_BEGIN, MARK_SCRIPTS_END, tables_md)

    if args.write:
        if content != new_content:
            # Backup
            if readme_path.exists():
                backup = readme_path.with_suffix(readme_path.suffix + ".bak")
                backup.write_text(content, encoding="utf-8")
            readme_path.write_text(new_content, encoding="utf-8")
            print(f"README updated → {readme_path}")
        else:
            print("README unchanged (no differences).")
    else:
        # Preview to stdout
        print("──────────────────── AUTO-OVERVIEW (preview) ────────────────────")
        print(overview_md)
        print("\n──────────────────── AUTO-SCRIPTS (preview) ─────────────────────")
        print(tables_md)

if __name__ == "__main__":
    main()
