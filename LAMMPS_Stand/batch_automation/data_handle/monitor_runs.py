#!/usr/bin/env python3
"""
Minimal Hall Monitor (with current stage via log tail)

Columns:
- run_id       : directory name (e.g., run.A.x2.r1)
- run_dir      : path relative to launch directory (PWD)
- state        : FINISHED (if any *.flag), ERROR (if any error flag), else RUNNING
- has_flag     : "1" if any MON_FLAG_GLOB matches in run dir, else "0"
- current_stage: last detected "Starting ..." stage from newest log (equilibration, production, ...)

Environment variables (override with qsub -v ...):
  MON_ROOT_DIR            default: $PWD/run_outputs
  MON_RUN_GLOB            default: run.*
  MON_FLAG_GLOB           default: *.flag
  MON_ERROR_FLAGS         default: run_error.flag,run_failed.flag
  MON_CSV_PATH            default: run_states.csv
  MON_INTERVAL_MIN        default: 60
  MON_STOP_WHEN_ALL_DONE  default: 1
  MON_MAX_HOURS           default: 240
  MON_INITIAL_DELAY_MIN   default: 0
  MON_LOG_GLOBS           default: job_output_*.txt,*.out,log.lammps,*.log
  MON_TAIL_BYTES          default: 1048576 (1 MiB) – how much of the log to scan from the end
"""

import os
import csv
import glob
import time
import signal
from datetime import datetime
from typing import List, Optional, Tuple

# ---------------- Config ----------------
ROOT_DIR            = os.environ.get("MON_ROOT_DIR", os.path.join(os.getcwd(), "run_outputs"))
RUN_GLOB            = os.environ.get("MON_RUN_GLOB", "run.*")
FLAG_GLOB           = os.environ.get("MON_FLAG_GLOB", "*.flag")
CSV_PATH            = os.environ.get("MON_CSV_PATH", "run_states.csv")
INTERVAL_MIN        = int(os.environ.get("MON_INTERVAL_MIN", "60"))
STOP_WHEN_ALL_DONE  = os.environ.get("MON_STOP_WHEN_ALL_DONE", "1").lower() in ("1","true","yes")
MAX_HOURS           = float(os.environ.get("MON_MAX_HOURS", "240"))
INITIAL_DELAY_MIN   = int(os.environ.get("MON_INITIAL_DELAY_MIN", "0"))

ERROR_FLAG_NAMES    = [x.strip() for x in os.environ.get("MON_ERROR_FLAGS", "run_error.flag,run_failed.flag").split(",") if x.strip()]
LOG_GLOBS           = [x.strip() for x in os.environ.get("MON_LOG_GLOBS", "job_output_*.txt,*.out,log.lammps,*.log").split(",") if x.strip()]
TAIL_BYTES          = int(os.environ.get("MON_TAIL_BYTES", str(1024 * 1024)))  # 1 MiB

# Map substrings in "Starting ..." lines to canonical stage names
STAGE_CANON_MAP = {
    "minimization": "minimization",
    "relaxation":   "relaxation",
    "nvt thaw":     "nvt_thaw",
    "thaw":         "nvt_thaw",
    "equil":        "equilibration",
    "equilibrat":   "equilibration",
    "production":   "production",
}

stop_now = False
def _handle_sigterm(signum, frame):
    global stop_now
    stop_now = True
signal.signal(signal.SIGTERM, _handle_sigterm)
signal.signal(signal.SIGINT, _handle_sigterm)

# ---------------- Helpers ----------------
def list_runs(root: str, pat: str) -> List[str]:
    """Absolute paths to run directories under root matching pat (recursive)."""
    # Search recursively under root; match terminal dirnames with pattern
    runs: List[str] = []
    for base, dirs, _files in os.walk(root):
        for d in dirs:
            if glob.fnmatch.fnmatch(d, pat):
                p = os.path.join(base, d)
                if os.path.isdir(p):
                    runs.append(p)
    runs.sort()
    return runs

def has_any_flag(run_dir_abs: str, pattern: str) -> bool:
    return any(glob.glob(os.path.join(run_dir_abs, pattern)))

def has_error_flag(run_dir_abs: str, names: List[str]) -> bool:
    for name in names:
        if os.path.isfile(os.path.join(run_dir_abs, name)):
            return True
    return False

def find_latest_log(run_dir: str) -> Optional[str]:
    """Find newest log file among LOG_GLOBS in run_dir; return path or None."""
    candidates: List[Tuple[float, str]] = []
    for pat in LOG_GLOBS:
        for p in glob.glob(os.path.join(run_dir, pat)):
            try:
                candidates.append((os.path.getmtime(p), p))
            except FileNotFoundError:
                pass
    if not candidates:
        return None
    candidates.sort(reverse=True)
    return candidates[0][1]

def tail_read(path: str, max_bytes: int) -> str:
    """Read up to last max_bytes from file (text, tolerate encoding issues)."""
    try:
        size = os.path.getsize(path)
        start = max(0, size - max_bytes)
        with open(path, "rb") as f:
            f.seek(start)
            data = f.read()
        return data.decode(errors="replace")
    except Exception:
        return ""

def detect_current_stage_from_text(text: str) -> str:
    """
    Look for the *last* occurrence of a line containing 'starting' and a known keyword.
    Return canonical stage or "".
    """
    current = ""
    for line in text.splitlines():
        s = line.strip().lower()
        if "starting" not in s:
            continue
        # find first matching keyword in priority order as written above
        for key, canon in STAGE_CANON_MAP.items():
            if key in s:
                current = canon
                break
    return current

def detect_current_stage(run_dir_abs: str) -> str:
    logp = find_latest_log(run_dir_abs)
    if not logp:
        return ""
    tail = tail_read(logp, TAIL_BYTES)
    return detect_current_stage_from_text(tail)

def scan_once():
    rows = []
    for run_dir_abs in list_runs(ROOT_DIR, RUN_GLOB):
        run_id = os.path.basename(run_dir_abs)
        run_dir_rel = os.path.relpath(run_dir_abs, os.getcwd())

        finished = has_any_flag(run_dir_abs, FLAG_GLOB)
        errored  = has_error_flag(run_dir_abs, ERROR_FLAG_NAMES)

        # If you want ERROR to trump FINISHED, swap the order below.
        state = "FINISHED" if finished else ("ERROR" if errored else "RUNNING")

        current_stage = detect_current_stage(run_dir_abs)

        rows.append({
            "run_id": run_id,
            "run_dir": run_dir_rel,
            "state": state,
            "has_flag": "1" if finished else "0",
            "current_stage": current_stage,
        })
    return rows

def write_csv(rows):
    os.makedirs(os.path.dirname(CSV_PATH) or ".", exist_ok=True)
    fieldnames = ["run_id", "run_dir", "state", "has_flag", "current_stage"]
    with open(CSV_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

# ---------------- Main ----------------
def main():
    if INITIAL_DELAY_MIN > 0:
        print(f"[monitor] Sleeping {INITIAL_DELAY_MIN} minute(s) before first scan...")
        for _ in range(INITIAL_DELAY_MIN * 60):
            if stop_now: return 0
            time.sleep(1)

    print(f"[run_monitor] root={ROOT_DIR} pattern={RUN_GLOB} flag_glob={FLAG_GLOB} interval={INTERVAL_MIN}min")
    print(f"[run_monitor] log_globs={','.join(LOG_GLOBS)} tail_bytes={TAIL_BYTES}")
    start = time.time()

    while True:
        rows = scan_once()
        write_csv(rows)

        total     = len(rows)
        finished  = sum(1 for r in rows if r["state"] == "FINISHED")
        errored   = sum(1 for r in rows if r["state"] == "ERROR")
        timestamp = datetime.now().isoformat(timespec="seconds")
        print(f"[{timestamp}] Scanned {total} runs -> FINISHED={finished}, ERROR={errored}. CSV: {CSV_PATH}", flush=True)

        if total > 0 and STOP_WHEN_ALL_DONE and finished == total:
            print("All runs finished. Exiting monitor.")
            return 0

        if (time.time() - start) / 3600.0 >= MAX_HOURS:
            print(f"Max hours ({MAX_HOURS}) reached. Exiting monitor.")
            return 0

        if stop_now:
            print("Received stop signal. Exiting monitor.")
            return 0

        # Sleep until next tick
        for _ in range(INTERVAL_MIN * 60):
            if stop_now:
                print("Received stop signal during sleep. Exiting monitor.")
                return 0
            time.sleep(1)

if __name__ == "__main__":
    raise SystemExit(main())
