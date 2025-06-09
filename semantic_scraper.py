#!/usr/bin/env python3
"""
semantic_scraper.py
────────────────────────────────────────────────────────────────────────
A beginner-friendly pipeline that

1. queries the **Semantic Scholar Graph API** for open-access (OA) papers,
2. downloads their PDFs (robustly and resumably),
3. extracts every table the PDF parser can detect,
4. keeps only tables that look numeric (→ e.g. diesel / elemental analysis),
5. writes a run summary + full log so you can audit failures.

The script prints **only tqdm progress bars** to the terminal; everything
else (warnings, retries, bad URLs, etc.) is written to *scrape.log* in the
output directory.

DEPENDENCIES  (conda-forge versions recommended)
─────────────────────────────────────────────────
  requests, pandas, tqdm, PyYAML
  camelot-py, ghostscript, opencv, tk         ← Camelot backend
  pdfplumber   (optional, OCR fallback)

Example install:
    conda install -c conda-forge requests pandas tqdm pyyaml \
                  camelot-py ghostscript opencv tk pdfplumber
────────────────────────────────────────────────────────────────────────

Author: Thomas J. Burton – Savoie Research Group, University of Notre Dame
Updated: 2025-06-06
License: MIT
────────────────────────────────────────────────────────────────────────
Usage:
python semantic_scraper.py -q 'Keyword Search" -n 999 -o output_dir -s tabs
"""

# ╭────────────────── USER-TUNABLE CONSTANTS ──────────────────────────╮
NUMERIC_RATIO_MIN = 0.20   # ≥ 20 % of data cells contain digits → keep
NUMERIC_CELLS_MIN = 9      # minimum absolute number of numeric cells
REQUEST_DELAY     = 3      # s between anonymous API calls
MAX_DL_THREADS    = 6      # parallel downloads; tweak to bandDownlwidth
# ╰────────────────────────────────────────────────────────────────────╯

import argparse
import html as ihtml
import logging
import os
import re
import sys
import textwrap
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import requests
from tqdm import tqdm

# ────────────────────────── Silence noisy warnings ───────────────────
warnings.filterwarnings(
    "ignore",
    message="CropBox missing from /Page, defaulting to MediaBox"
)

# ────────────────────────── Required / optional libs ─────────────────
try:
    import camelot  # table extractor for vector PDFs
except ImportError:
    sys.exit(
        "❌  Camelot missing.\n"
        "    Install: conda install -c conda-forge camelot-py ghostscript opencv tk"
    )

try:
    import pdfplumber  # OCR backup for image-only tables
    PDFPLUMBER_OK = True
except ImportError:
    PDFPLUMBER_OK = False

# ────────────────────────── API constants ────────────────────────────
S2_API_ROOT   = "https://api.semanticscholar.org/graph/v1"
DEFAULT_FIELDS = (
    "title,year,authors,url,openAccessPdf,abstract,publicationTypes"
)

HEADERS_API = (
    {"x-api-key": os.getenv("S2_API_KEY")}
    if os.getenv("S2_API_KEY")
    else {}
)

# regex used in numeric-heuristic
DIGIT_TOKEN = re.compile(r"(\d|\bppm\b|\bwt\b|\bvol\b)", re.I)

# ╭──────────────────── Helper: API GET with 429 handling ─────────────╮
def polite_get(url: str, params: Dict[str, str], max_retries: int = 5) -> Dict:
    """
    Thin wrapper around `requests.get` that
    • sleeps REQUEST_DELAY seconds between calls,
    • backs off exponentially on HTTP 429,
    • logs every retry to scrape.log,
    • raises for any other HTTP error.
    """
    delay = REQUEST_DELAY
    for attempt in range(1, max_retries + 1):
        time.sleep(delay)
        resp = requests.get(url, headers=HEADERS_API, params=params, timeout=30)

        if resp.status_code == 200:                 # happy path
            return resp.json()

        if resp.status_code == 429:                 # rate-limited
            wait = int(resp.headers.get("Retry-After", "60"))
            logging.warning(f"429 rate-limit: sleeping {wait}s  (retry {attempt})")
            print("\n❗️  Rate limit hit, sleeping for", wait, "seconds...\n")
            time.sleep(wait)
            delay *= 2
            continue

        resp.raise_for_status()                     # other 4xx / 5xx

    raise RuntimeError("Exceeded retry budget for API 429s")

# ╭────────────────────── Search helper (bulk + paged) ───────────────╮
def search_papers(query: str, limit: int) -> List[Dict]:
    """
    Return up to `limit` paper JSON blobs *that have an OA PDF*.
    Fast path  : `/paper/search/bulk` when limit ≤ 1000.
    Slow path  : classic pagination (offset), stops at offset 1000.
    """
    collected: List[Dict] = []

    # ---------- bulk endpoint ----------
    if limit <= 1000:
        payload = {"query": query, "limit": limit, "fields": DEFAULT_FIELDS}
        logging.info(f"Bulk search  query={query!r}")
        for p in polite_get(f"{S2_API_ROOT}/paper/search/bulk", payload)["data"]:
            if p.get("openAccessPdf", {}).get("url"):
                collected.append(p)
        return collected

    # ---------- paged search ----------
    offset, page = 0, 100
    while len(collected) < limit and offset < 1000:
        payload = {
            "query": query,
            "limit": min(page, limit - len(collected)),
            "offset": offset,
            "fields": DEFAULT_FIELDS,
        }
        if offset == 0:
            logging.info(f"Paged search URL params: {payload}")
        try:
            page_data = polite_get(f"{S2_API_ROOT}/paper/search", payload)["data"]
        except requests.HTTPError as exc:
            logging.warning(f"Stopping search at offset {offset}: {exc}")
            break

        for p in page_data:
            if p.get("openAccessPdf", {}).get("url"):
                collected.append(p)
                if len(collected) == limit:
                    break
        offset += page
        if not page_data:
            break

    return collected

# ╭──────────────────── PDF utilities ───────────────────────────────╮
def is_valid_pdf(path: Path, *, min_bytes: int = 5_000) -> bool:
    """Cheap check that `path` points to a real PDF (EOF marker present)."""
    if not path.exists() or path.stat().st_size < min_bytes:
        return False
    with path.open("rb") as fh:
        try:
            fh.seek(-32, os.SEEK_END)
        except OSError:
            return False
        tail = fh.read()
        return b"%%EOF" in tail or b"endobj" in tail

def sanitize(name: str) -> str:
    """File-system-safe slug of a title (max 25 chars)."""
    return re.sub(r'[\\/:*?"<>|]+', "_", name)[:25]

# ╭──────────────────── Robust PDF downloader ───────────────────────╮
def download_pdf(
    url: str,
    dest: Path,
    referer: str,
    max_retries: int = 3,
) -> Tuple[bool, str]:
    """
    Download a PDF (with landing-page & MDPI hacks).

    Returns
    -------
    (True, "")              → success or already on disk
    (False, "reason text")  → all retries failed
    """
    if is_valid_pdf(dest):
        return True, ""

    hdr = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/125.0.0.0 Safari/537.36"
        ),
        "Referer": referer,
    }
    err_msg = ""
    for attempt in range(1, max_retries + 1):
        try:
            with requests.get(
                url,
                headers=hdr,
                stream=True,
                timeout=60,
                allow_redirects=True,
            ) as resp:

                # ‣ Landing page returned instead of PDF
                ctype = resp.headers.get("Content-Type", "").lower()
                if ctype.startswith("text/html"):
                    html = resp.text
                    # meta-refresh
                    m = re.search(
                        r'http-equiv=["\']refresh["\'].*?url=(.*?)["\']',
                        html,
                        flags=re.I,
                    )
                    if m:
                        url = ihtml.unescape(m.group(1))
                        logging.info(f"meta-refresh → {url}")
                        time.sleep(1)
                        continue
                    # direct <a> to pdf
                    a = re.search(r'href=["\'](.*?\.pdf.*?)["\']', html, flags=re.I)
                    if a:
                        href = ihtml.unescape(a.group(1))
                        url = requests.compat.urljoin(resp.url, href)
                        logging.info(f"HTML link → {url}")
                        time.sleep(1)
                        continue
                    raise requests.HTTPError("HTML page, no PDF link")

                # ‣ MDPI occasionally 403s without ?download=1
                if resp.status_code == 403 and "mdpi.com" in url:
                    alt = url.split("?")[0] + "?download=1"
                    logging.info(f"MDPI 403 → retry {alt}")
                    url = alt
                    time.sleep(1)
                    continue

                # ‣ Save the stream to disk
                resp.raise_for_status()
                with dest.open("wb") as fh:
                    for chunk in resp.iter_content(chunk_size=2 << 14):
                        fh.write(chunk)

                if is_valid_pdf(dest):
                    return True, ""
                raise requests.HTTPError("File saved but invalid PDF")

        except requests.RequestException as exc:
            err_msg = f"{type(exc).__name__}: {exc}"
            logging.warning(f"({attempt}/{max_retries}) {err_msg}")
            time.sleep(2 * attempt)

    return False, err_msg

# ╭────────────────── Heuristic to accept / reject tables ───────────╮
def table_is_useful(df: pd.DataFrame) -> bool:
    """Reject one-column or non-numeric tables."""
    if df.shape[1] < 2:
        return False
    data_rows = df.iloc[1:] if len(df) > 1 else df
    total, numericish = 0, 0
    for cell in data_rows.values.flatten():
        s = str(cell)
        if len(s) > 20 or s.startswith("("):  # long footnote / DOI
            continue
        total += 1
        if DIGIT_TOKEN.search(s):
            numericish += 1
    if total == 0:
        return False
    ratio = numericish / total
    return ratio >= NUMERIC_RATIO_MIN and numericish >= NUMERIC_CELLS_MIN

# ╭──────────────────── Worker: extract tables per PDF ───────────────╮
def extract_tables_worker(
    args: Tuple[Path, int, Path, str, str]
) -> Tuple[int, int]:
    """
    Extract every table and write as CSV with metadata header lines.

    Parameters
    ----------
    pdf_path     : Path
    idx          : int   Paper index (1-based)
    tables_dir   : Path  Common output folder
    paper_title  : str   Human-readable title
    paper_url    : str   Semantic Scholar canonical URL

    Returns
    -------
    (idx, raw_tables_saved)
    """
    pdf_path, idx, tables_dir, paper_title, paper_url = args
    import camelot, pandas as pd  # local import inside subprocess

    saved = 0

    def _save(df: pd.DataFrame) -> None:
        nonlocal saved
        out_file = tables_dir / f"{idx:03d}_T{saved+1:02d}.csv"
        with out_file.open("w", encoding="utf-8", newline="") as fh:
            fh.write(f"# title: {paper_title}\n")
            fh.write(f"# url:   {paper_url}\n")
            df.to_csv(fh, index=False)
        saved += 1

    # 1) Camelot (vector PDFs)
    for flavor in ("lattice", "stream"):
        try:
            for tbl in camelot.read_pdf(str(pdf_path), flavor=flavor, pages="all"):
                _save(tbl.df)
        except Exception:
            pass

    # 2) pdfplumber fallback (image PDFs)
    if saved == 0 and PDFPLUMBER_OK:
        try:
            import pdfplumber
            with pdfplumber.open(str(pdf_path)) as pdf:
                for page in pdf.pages:
                    for raw in page.extract_tables():
                        _save(pd.DataFrame(raw))
        except Exception:
            pass

    return idx, saved

# ╭───────────────────────────── main() ──────────────────────────────╮
def main() -> None:
    # ── CLI parsing ──────────────────────────────────────────────────
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Scrape numeric composition tables from OA PDFs.",
    )
    ap.add_argument(
        "-n",
        "--num",
        type=int,
        default=999,
        help="Number of papers to process (≤1000 per query)",
    )
    ap.add_argument(
        "-q",
        "--query",
        default="Rice Pyrolysis",
        help="Search string (use quotes for multi-word)",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        default="rice_pyrolysis_tables",
        help="Output directory (can be on external drive)",
    )
    args = ap.parse_args()

    # ── Prepare folders ──────────────────────────────────────────────
    root       = Path(args.outdir).resolve()
    tables_dir = root / "tables"
    junk_dir   = tables_dir / "_junk"
    for d in (root, tables_dir, junk_dir):
        d.mkdir(parents=True, exist_ok=True)

    # ── File-based logging ───────────────────────────────────────────
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)s  %(message)s",
        handlers=[logging.FileHandler(root / "scrape.log", encoding="utf-8")],
    )
    logging.info(f"RUN  query={args.query!r}  target={args.num}")

    # ── 1. Metadata search ───────────────────────────────────────────
    print(f"🔎  Query={args.query!r}  Target={args.num}")
    papers = search_papers(args.query, args.num)
    print(f"   ↳ OA hits: {len(papers)}\n")
    logging.info(f"OA papers returned: {len(papers)}")

    # Build job lists / summary skeleton
    dl_jobs: List[Tuple[int, Dict, Path]] = []
    summary_rows: List[Dict] = []
    for idx, paper in enumerate(papers, 1):
        pdir = root / f"{idx:03d}_{sanitize(paper.get('title','untitled'))}"
        pdir.mkdir(exist_ok=True)
        dl_jobs.append((idx, paper, pdir))

    # ── 2. Parallel PDF downloads (ThreadPool) ───────────────────────
    def thread_dl(job):
        idx, paper, pdir = job
        pdf_url  = paper["openAccessPdf"]["url"]
        pdf_path = pdir / "paper.pdf"
        ok, err  = download_pdf(pdf_url, pdf_path, paper["url"])
        return idx, ok, err, pdf_path

    extract_jobs: List[Tuple[Path, int, Path, str, str]] = []

    with ThreadPoolExecutor(max_workers=MAX_DL_THREADS) as pool:
        for idx, ok, err, pdf_path in tqdm(
            pool.map(thread_dl, dl_jobs),
            total=len(dl_jobs),
            ncols=80,
            desc="Downloading",
        ):
            paper = papers[idx - 1]
            summary_rows.append(
                {
                    "title": paper["title"],
                    "url": paper["url"],
                    "pdf_status": "ok" if ok else "failed",
                    "error": err,
                    "tables_raw": 0,
                    "tables_kept": 0,
                }
            )
            if ok:
                # pass title & URL to extraction worker
                extract_jobs.append(
                    (pdf_path, idx, tables_dir, paper["title"], paper["url"])
                )
            else:
                logging.warning(f"PDF FAILED  {paper['url']}  » {err}")

    # ── 3. Table extraction (ProcessPool) ────────────────────────────
    if extract_jobs:
        with ProcessPoolExecutor(max_workers=os.cpu_count()) as pool:
            for idx, raw in tqdm(
                pool.map(extract_tables_worker, extract_jobs),
                total=len(extract_jobs),
                ncols=80,
                desc="Extracting",
            ):
                summary_rows[idx - 1]["tables_raw"] = raw

    # ── 4. Junk filter pass ─────────────────────────────────────────
    kept_by_idx: Dict[int, int] = {}
    for csv in tables_dir.glob("*.csv"):
        # skip comment header lines when re-reading
        df = pd.read_csv(csv, dtype=str, keep_default_na=False, comment="#")
        i  = int(csv.stem.split("_")[0])
        if table_is_useful(df):
            kept_by_idx[i] = kept_by_idx.get(i, 0) + 1
        else:
            csv.rename(junk_dir / csv.name)

    for i, row in enumerate(summary_rows, 1):
        row["tables_kept"] = kept_by_idx.get(i, 0)

    logging.info(
        f"Tables kept: {sum(kept_by_idx.values())}  "
        f"Junk moved: {len(list(junk_dir.glob('*.csv')))}"
    )

    # ── 5. Write summary sheets ─────────────────────────────────────
    pd.DataFrame(summary_rows).to_csv(root / "summary.csv", index=False)
    pd.DataFrame(summary_rows).to_excel(root / "summary.xlsx", index=False)

    # ── 6. Remove empty/failed paper folders ────────────────────────
    import shutil

    removed = 0
    for idx, row in enumerate(summary_rows, 1):
        if row["pdf_status"] != "failed":
            continue
        pdir = root / f"{idx:03d}_{sanitize(row['title'])}"
        try:
            shutil.rmtree(pdir)
            removed += 1
        except (FileNotFoundError, OSError):
            pass
    logging.info(f"Empty dirs removed: {removed}")

    # ── 7. Final console / log messages ─────────────────────────────
    print(f"\n✅  DONE — summary & log in {root}")
    if not PDFPLUMBER_OK:
        logging.info("pdfplumber not installed (OCR fallback skipped)")
    logging.info("RUN finished\n")

# ──────────────────────────────── run ────────────────────────────────
if __name__ == "__main__":
    main()
