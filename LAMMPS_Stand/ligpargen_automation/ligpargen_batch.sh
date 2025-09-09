#!/usr/bin/env bash
set -euo pipefail

XLSX="Electrolytes.xlsx"
SCRAPER="scrape_py.py"
OUTDIR="output_lmp"
MAPCSV="ligpargen_mapping.csv"

# ---------- sanity checks ----------
[[ -f "$XLSX" ]]    || { echo "[ERROR] Missing $XLSX (run the XLSX-maker first)"; exit 1; }
[[ -f "$SCRAPER" ]] || { echo "[ERROR] Missing $SCRAPER"; exit 1; }

# ---------- ensure chromedriver is available as ./chromedriver ----------
need_link=true
if [[ -x "./chromedriver" ]]; then
  echo "[INFO] Found ./chromedriver"
  need_link=false
fi

if $need_link; then
  echo "[INFO] Locating chromedriver on this system..."
  cdpath=""
  # try common locations / PATH / Homebrew
  if command -v chromedriver >/dev/null 2>&1; then
    cdpath="$(command -v chromedriver)"
  elif [[ -x "/usr/local/bin/chromedriver" ]]; then
    cdpath="/usr/local/bin/chromedriver"
  elif [[ -x "/opt/homebrew/bin/chromedriver" ]]; then
    cdpath="/opt/homebrew/bin/chromedriver"
  else
    # last try: brew --prefix (macOS Homebrew)
    if command -v brew >/dev/null 2>&1; then
      bp="$(brew --prefix 2>/dev/null || true)"
      if [[ -n "${bp:-}" && -x "$bp/bin/chromedriver" ]]; then
        cdpath="$bp/bin/chromedriver"
      fi
    fi
  fi

  if [[ -n "${cdpath:-}" && -x "$cdpath" ]]; then
    echo "[INFO] Using chromedriver at: $cdpath"
    ln -sfn "$cdpath" "./chromedriver"
    chmod +x "./chromedriver" || true
  else
    echo "[ERROR] Could not find a chromedriver binary."
    echo "        Options to fix:"
    echo "          • macOS (Homebrew):   brew install --cask chromedriver"
    echo "          • Linux (APT):        sudo apt-get install chromium-chromedriver"
    echo "          • Or download a matching version from: https://chromedriver.chromium.org/"
    echo "        After install, re-run this script."
    exit 1
  fi
fi

# ---------- run the scraper (UNCHANGED code expects ./chromedriver) ----------
echo "[INFO] Running scraper: $SCRAPER"
python3 "$SCRAPER"

# ---------- build mapping CSV from the XLSX ----------
echo "[INFO] Building mapping CSV: $MAPCSV"
python3 - <<'PYEOF'
import os, pandas as pd

xlsx = "Electrolytes.xlsx"
outdir = "output_lmp"
mapcsv = "ligpargen_mapping.csv"

df = pd.read_excel(xlsx)
if "rdkit_smiles" not in df.columns:
    raise SystemExit("[ERROR] 'rdkit_smiles' column not found in Electrolytes.xlsx")

def lmp_name(smiles: str) -> str:
    # Matches scrape_py.py behavior: smiles.replace('/', '_') + ".lmp"
    return f"{str(smiles).replace('/', '_')}.lmp"

df["smiles"]   = df["rdkit_smiles"].astype(str)
df["lmp_file"] = df["smiles"].map(lmp_name)
df["exists"]   = df["lmp_file"].apply(lambda f: os.path.isfile(os.path.join(outdir, f)))

df[["smiles","lmp_file","exists"]].to_csv(mapcsv, index=False)
ok = int(df["exists"].sum())
print(f"[INFO] Wrote {len(df)} rows -> {mapcsv}  (found {ok} files)")
PYEOF

echo "[INFO] Done."
echo " - LAMMPS files:   $OUTDIR/"
echo " - Mapping CSV:    $MAPCSV"
