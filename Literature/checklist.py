#!/usr/bin/env python3
"""
due_today.py
------------
Check which papers from reading_log.csv are due for re-reading.

Default: shows papers read exactly 14 days ago and not yet revisited.
Optional: use --overdue to also list all older entries not yet revisited.

Accepted date formats include:
  - 2025-11-03
  - 03NOV2025
  - 03-Nov-2025
  - 03/11/2025
"""

import pandas as pd
import datetime as dt
import re
import argparse


def parse_date_any(s):
    """Try multiple date formats (case-insensitive)."""
    if pd.isna(s):
        return pd.NaT
    s = str(s).strip()

    ts = pd.to_datetime(s, errors="coerce")
    if not pd.isna(ts):
        return ts.normalize()

    # DDMMMYYYY (e.g., 03NOV2025)
    m = re.fullmatch(r"(\d{2})([A-Za-z]{3})(\d{4})", s)
    if m:
        day, mon3, year = m.groups()
        ts = pd.to_datetime(f"{day}{mon3.title()}{year}", format="%d%b%Y", errors="coerce")
        if not pd.isna(ts):
            return ts.normalize()

    # Common explicit formats
    for fmt in ("%d-%b-%Y", "%d/%m/%Y", "%m/%d/%Y", "%Y-%m-%d"):
        ts = pd.to_datetime(s, format=fmt, errors="coerce")
        if not pd.isna(ts):
            return ts.normalize()

    return pd.NaT


def main():
    p = argparse.ArgumentParser(description="Check for papers due for re-reading.")
    p.add_argument("--days", type=int, default=14, help="Days since first read (default=14)")
    p.add_argument("--overdue", action="store_false",default=True, help="Also show older, unrevisited papers")
    args = p.parse_args()

    df = pd.read_csv("reading_log.csv")
    if "Date" not in df.columns:
        raise SystemExit("reading_log.csv must have a 'Date' column")

    df["_Date_parsed"] = df["Date"].apply(parse_date_any)
    if df["_Date_parsed"].isna().all():
        raise SystemExit("Could not parse any dates in 'Date' column. Check formats.")

    today = pd.Timestamp(dt.date.today())
    target_date = (today - pd.Timedelta(days=args.days)).date()

    revisited_col = "Revisited" if "Revisited" in df.columns else None
    if revisited_col:
        revisited_true = df[revisited_col].astype(str).str.lower().isin(["yes", "true", "1", "y"])
        df = df[~revisited_true]

    df["_Date_only"] = df["_Date_parsed"].dt.date

    if args.overdue:
        due = df[df["_Date_only"] <= target_date]
    else:
        due = df[df["_Date_only"] == target_date]

    print(f"📅 Papers {'overdue or due' if args.overdue else 'due'} (read on/before {target_date}):")
    if due.empty:
        print("None found")
        return

    cols = [c for c in ["Date", "Paper", "Relevance", "Revisited"] if c in due.columns]
    print(due[cols].to_string(index=False))


if __name__ == "__main__":
    main()
