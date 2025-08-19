#!/usr/bin/env python3
"""
update_readme.py
───────────────────────────────────────────────────────────────────────
Auto-refreshes README sections by scanning the repository for scripts.

Features
• Recursively discovers .py/.sh scripts, skipping junk (e.g., .git, venvs)
• Extracts:
  – Title (filename), relative path
  – Language (Python/Bash)
  – Description (module docstring or header comment)
  – Usage (from docstring, header, or optional --help capture)
  – Guessed dependencies (imports / external commands)
  – Lines of code, last modified time
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
  python update_readme.py --run-help    # try script -h/--help capture
  python update_readme.py --write-index # emit script_index.json/.csv

Recommended exclusions for speed/sanity:
  --exclude-dirs .git .venv build dist __pycache__ node_modules .mypy_cache
"""

from __future__ import annotations

import argparse
import ast
import csv
import datetime as dt
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List, Dict, Optional, Tuple

# ─────────────────────────────────────────────────────────────────────
# Config (defaults can be overridden via CLI)
# ─────────────────────────────────────────────────────────────────────
READ_ME = Path("README.md")
INCLUDE_EXT = {".py", ".sh"}
DEFAULT_EXCLUDE_DIRS = {
    ".git", ".svn", ".hg", ".mypy_cache", "__pycache__",
    ".venv", "venv", "env", "build", "dist", ".ipynb_checkpoints",
    ".pytest_cache", ".tox", ".idea", ".vscode", ".DS_Store",
}
IGNORE_HIDDEN = True
HELP_TIMEOUT = 7  # seconds per script when --run-help is set

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

PY_STDLIB_APPROX = {
    # Minimal set for filtering. This is intentionally small; we only
    # elevate *non*-stdlib as "guessed deps".
    "os","sys","re","json","csv","pathlib","argparse","textwrap","typing",
    "subprocess","time","datetime","math","itertools","logging","collections",
    "functools","shutil","ast","tempfile","warnings","dataclasses","statistics",
    "html","urllib","enum","hashlib","gzip","bz2","io","glob","gc","pickle",
}

BASH_COMMON_CMDS = {
    "ffmpeg","awk","sed","grep","cut","tr","sort","uniq","head","tail","xargs",
    "jq","curl","wget","tar","gzip","rg","fd","perl","python","python3","bash",
    "zsh","sh","tee","mktemp","convert","identify","parallel","gs"
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
    guessed_deps: List[str]
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
            # ensure excluded top-level fragments aren’t inside the path
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

def first_paragraph(s: str, max_chars: int = 500) -> str:
    s = s.strip().splitlines()
    block = []
    for line in s:
        if line.strip().startswith(("Usage", "USAGE")):
            break
        if line.strip().startswith("---"):
            break
        block.append(line)
    out = " ".join(line.strip() for line in block if line.strip())
    return (out[:max_chars] + "…") if len(out) > max_chars else out

def extract_usage_block(text: str) -> str:
    # Very simple heuristic: look for a "Usage" section or lines with "python <file>".
    usage = []
    m = re.search(r"(?i)usage\s*:?(.*)", text)
    if m:
        # capture from "Usage" to the next blank line / dashed rule
        start = m.start()
        tail = text[start:].split("\n\n", 1)[0]
        # keep short, de-indent
        for line in tail.splitlines():
            line = line.strip().strip("`")
            if not line:
                break
            usage.append(line)
    else:
        # fallback: grep for obvious CLI pattern
        for line in text.splitlines():
            if re.search(r"\bpython3?\s+\S+\.py\b", line) or re.search(r"\b\.\/\S+\.sh\b", line):
                usage.append(line.strip())
    usage = [u for u in usage if len(u) < 160][:3]
    return "\n".join(usage)

def parse_python_info(path: Path, text: str, run_help: bool) -> Tuple[str, str, List[str]]:
    """
    Returns (description, usage, guessed_deps)
    """
    desc = ""
    usage = ""
    deps: set[str] = set()

    try:
        mod = ast.parse(text)
        doc = ast.get_docstring(mod) or ""
        desc = first_paragraph(doc) or "Python script"
        usage = extract_usage_block(doc)
        for node in ast.walk(mod):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                name = node.names[0].name.split(".")[0]
                if name not in PY_STDLIB_APPROX:
                    deps.add(name)
    except Exception:
        # Fallback to header comments
        desc = first_paragraph(text) or "Python script"

    if run_help:
        help_text = try_help(path)
        if help_text:
            # Prefer explicit help usage if present
            usage = extract_usage_block(help_text) or usage
            # mine deps hints
            m = re.findall(r"(?i)requires?:\s*(.+)", help_text)
            for g in m:
                for token in re.split(r"[,\s]+", g):
                    t = token.strip()
                    if t and t not in PY_STDLIB_APPROX:
                        deps.add(t)

    # Curate known libs so we don’t list huge stdlib guesses
    curated = []
    for d in sorted(deps):
        # light filter for popular libs we actually use
        if d.lower() in {
            "rdkit","requests","pandas","tqdm","camelot","pdfplumber","xlsxwriter",
            "openpyxl","opencv","cv2","ase","matplotlib","pillow","pil","pyyaml",
        }:
            curated.append(d)
    return desc, usage, curated

def parse_bash_info(path: Path, text: str, run_help: bool) -> Tuple[str, str, List[str]]:
    """
    Returns (description, usage, guessed_deps)
    """
    # Header comment block (ignoring shebang)
    header = []
    lines = text.splitlines()
    i = 0
    if lines and lines[0].startswith("#!"):
        i = 1
    for line in lines[i:]:
        if line.strip().startswith("#"):
            header.append(line.strip().lstrip("#").strip())
        elif not line.strip():
            # allow a single blank in header
            header.append("")
        else:
            break
    header_text = "\n".join(header).strip()
    desc = first_paragraph(header_text) or "Bash script"
    usage = extract_usage_block(header_text)

    # Guess deps by scanning tokens that look like external commands
    deps = set()
    tokens = re.findall(r"\b[a-zA-Z0-9_\-\.]+\b", text)
    for t in tokens:
        if t in BASH_COMMON_CMDS:
            deps.add(t)

    if run_help:
        help_text = try_help(path)
        if help_text:
            usage = extract_usage_block(help_text) or usage
            # Also look for dependencies lines
            for m in re.findall(r"(?i)requires?:\s*(.+)", help_text):
                for token in re.split(r"[,\s]+", m):
                    token = token.strip()
                    if token:
                        deps.add(token)

    return desc, usage, sorted(deps)

def try_help(path: Path) -> str:
    """
    Attempt to capture --help/-h.
    NEVER run arbitrary scripts without the user's consent: gated by --run-help.
    """
    try:
        if path.suffix == ".py":
            cmd = [sys.executable, str(path), "-h"]
        else:
            cmd = [str(path), "--help"]
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=HELP_TIMEOUT, check=False, text=True
        )
        return proc.stdout.strip()
    except Exception:
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
        out.append("| Script | Lang | Description | Usage | Deps | LOC | Modified |")
        out.append("|---|---|---|---|---|---|---|")
        for s in sorted(by_cat[cat], key=lambda x: x.name.lower()):
            link = f"[`{s.name}`]({s.relpath})"
            desc = s.description.replace("|", "\\|")
            usage = s.usage.replace("|", "\\|") if s.usage else ""
            deps = ", ".join(s.guessed_deps)
            out.append(
                f"| {link} | {s.language} | {desc} | {usage} | {deps} | {s.loc} | {s.modified} |"
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
    ap = argparse.ArgumentParser(description="Auto-update README tables from repository scripts.")
    ap.add_argument("--repo-root", type=Path, default=Path("."), help="Repository root (default: .)")
    ap.add_argument("--readme", type=Path, default=READ_ME, help="README path (default: README.md)")
    ap.add_argument("--include-ext", nargs="*", default=list(INCLUDE_EXT), help="File extensions to include")
    ap.add_argument("--exclude-dirs", nargs="*", default=list(DEFAULT_EXCLUDE_DIRS), help="Directory names to exclude")
    ap.add_argument("--no-ignore-hidden", action="store_true", help="Do not skip hidden files/folders")
    ap.add_argument("--run-help", action="store_true", help="Try to capture --help/-h output (may execute scripts)")
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
        text = read_text(path)
        lang = "Python" if path.suffix == ".py" else "Bash"

        if lang == "Python":
            desc, usage, deps = parse_python_info(path, text, args.run_help)
        else:
            desc, usage, deps = parse_bash_info(path, text, args.run_help)

        info = ScriptInfo(
            path=path.relative_to(root),
            language=lang,
            description=desc or (f"{lang} script"),
            usage=usage,
            guessed_deps=deps,
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
