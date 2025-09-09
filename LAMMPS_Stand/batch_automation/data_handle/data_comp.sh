#!/usr/bin/env bash
set -euo pipefail

# ==============================
# Recursive Run Data Summarizer
# ==============================

# Root that contains nested run.* dirs (override with: RUN_ROOT=/path ./collect_run_data.sh)
RUN_ROOT="${RUN_ROOT:-run_outputs}"

# Output files
OUT_CSV="data_summary.csv"
OUT_XLSX="data_summary.xlsx"   # (not created here, but left for compatibility)

# CSV header
echo "run_dir,run_finished,job_started,job_ended,wall_time_s,hostname,uge_job_id,run_name,D_li_cm2_s,D_cl_cm2_s,my_D_li_cm2_s,my_D_cl_cm2_s,mda_D_li_cm2_s,mda_D_cl_cm2_s" > "$OUT_CSV"

# Simple CSV escaping for text fields
csvq() {
  local s="${1:-}"
  s="${s//\"/\"\"}"
  printf '"%s"' "$s"
}

# Find all run.* directories recursively under RUN_ROOT (relative paths from PWD)
mapfile -t RUN_DIRS < <(find "$RUN_ROOT" -type d -name 'run.*' -print0 | xargs -0 -I{} realpath --relative-to="$PWD" "{}" | sort)

# If no runs found, exit gracefully
if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
  echo "[info] No run.* directories found under '$RUN_ROOT'."
  echo "[ok] Summary written to $OUT_CSV"
  exit 0
fi

for dir in "${RUN_DIRS[@]}"; do
  [[ -d "$dir" ]] || continue

  flag="$dir/run_finished.flag"

  # Defaults
  run_finished="false"
  job_started=""
  job_ended=""
  wall_time_s=""
  hostname=""
  uge_job_id=""
  run_name=""
  # “current info” D values (from flag file, as in your original script)
  D_li=""
  D_cl=""
  # New: values parsed from diffusion_results.txt (newest in the run)
  my_D_li=""
  my_D_cl=""
  mda_D_li=""
  mda_D_cl=""

  # --------- Parse run_finished.flag (if present) ----------
  if [[ -f "$flag" ]]; then
    run_finished="true"

    job_started="$(sed -n 's/^Job started at:[[:space:]]*//p' "$flag" | head -n1)"
    job_ended="$(sed -n 's/^Job ended at:[[:space:]]*//p' "$flag" | head -n1)"
    wall_time_s="$(sed -n 's/^Total wall time:[[:space:]]*\([0-9][0-9]*\)[[:space:]]*seconds.*/\1/p' "$flag" | head -n1)"
    hostname="$(sed -n 's/^Hostname:[[:space:]]*//p' "$flag" | head -n1)"
    uge_job_id="$(sed -n 's/^UGE Job ID:[[:space:]]*//p' "$flag" | head -n1)"
    run_name="$(sed -n 's/^Run Name:[[:space:]]*//p' "$flag" | head -n1)"

    # “Current info” diffusion constants captured in the flag file (tolerant to extra text)
    D_li="$(sed -nE 's/^msd_li\.txt:[[:space:]]*D[[:space:]]*=[[:space:]]*([0-9.eE+-]+).*/\1/p' "$flag" | head -n1)"
    D_cl="$(sed -nE 's/^msd_cl\.txt:[[:space:]]*D[[:space:]]*=[[:space:]]*([0-9.eE+-]+).*/\1/p' "$flag" | head -n1)"
  fi

  # --------- Find the newest diffusion_results.txt under this run ----------
  # Use -print0 to handle spaces; sort by mtime descending, pick first
  mapfile -d '' -t DIFF_FILES < <(find "$dir" -type f -name 'diffusion_results.txt' -print0)
  if [[ ${#DIFF_FILES[@]} -gt 0 ]]; then
    # Get the newest by mtime
    newest_diff="$(printf '%s\0' "${DIFF_FILES[@]}" | xargs -0 ls -1t | head -n1)"

    # ---- Parse “my diffusion” (top lines like: msd_li.txt: D = 1.1539e-06 cm²/s) ----
    my_D_li="$(sed -nE 's/^msd_li\.txt:[[:space:]]*D[[:space:]]*=[[:space:]]*([0-9.eE+-]+).*/\1/p' "$newest_diff" | head -n1)"
    my_D_cl="$(sed -nE 's/^msd_cl\.txt:[[:space:]]*D[[:space:]]*=[[:space:]]*([0-9.eE+-]+).*/\1/p' "$newest_diff" | head -n1)"

    # ---- Parse MDAnalysis “[Result] X: ...  D = 1.212e-06 cm^2/s” lines ----
    # Accept both "cm²/s" and "cm^2/s"
    mda_D_li="$(sed -nE 's/^\[Result\][[:space:]]*Li:.*D[[:space:]]*=[[:space:]]*([0-9.eE+-]+)[[:space:]]*cm(\^?2|²)\/s.*/\1/p' "$newest_diff" | head -n1)"
    mda_D_cl="$(sed -nE 's/^\[Result\][[:space:]]*Cl:.*D[[:space:]]*=[[:space:]]*([0-9.eE+-]+)[[:space:]]*cm(\^?2|²)\/s.*/\1/p' "$newest_diff" | head -n1)"
  fi

  # --------- Write CSV row (one row per run) ----------
  {
    csvq "$dir";             echo -n ","
    csvq "$run_finished";    echo -n ","
    csvq "$job_started";     echo -n ","
    csvq "$job_ended";       echo -n ","
    csvq "$wall_time_s";     echo -n ","
    csvq "$hostname";        echo -n ","
    csvq "$uge_job_id";      echo -n ","
    csvq "$run_name";        echo -n ","
    # Numeric-like fields left unquoted (blank if not found)
    echo -n "${D_li},${D_cl},${my_D_li},${my_D_cl},${mda_D_li},${mda_D_cl}"
    echo
  } >> "$OUT_CSV"

done

echo "[ok] Summary written to $OUT_CSV"
[[ -f "$OUT_XLSX" ]] && echo "[ok] Excel written to $OUT_XLSX"
