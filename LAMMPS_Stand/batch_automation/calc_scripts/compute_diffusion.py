import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import sys
import os
import argparse

import argparse

def analyze_msd(filename, outdir=None, timestep_fs=0.5, fit_start=1000, fit_end=7500):
    try:
        # --- Load MSD data ---
        data = np.loadtxt(filename)
        time_fs = data[:, 0] * timestep_fs
        time_ps = time_fs / 1000.0  # ps
        msd = data[:, 1]

        # --- Linear fit ---
        slope, intercept, r, p, stderr = linregress(time_ps[fit_start:fit_end],
                                                   msd[fit_start:fit_end])
        D_cm2s = (slope * 1e-4) / 6.0  # Correct conversion

        print(f"{filename}: D = {D_cm2s:.4e} cm²/s  (slope={slope})")

        # --- Plot ---
        plt.figure(figsize=(6, 4))
        plt.plot(time_ps, msd, label="MSD", lw=2)
        fit_line = intercept + slope * time_ps[fit_start:fit_end]
        plt.plot(time_ps[fit_start:fit_end], fit_line, "r--", label="Linear fit")
        plt.xlabel("Time (ps)")
        plt.ylabel("MSD (Å²)")
        plt.title(f"{filename}\nD = {D_cm2s:.2e} cm²/s")
        plt.legend()
        plt.tight_layout()

        # --- Output file ---
        base = os.path.splitext(os.path.basename(filename))[0]
        if outdir:
            os.makedirs(outdir, exist_ok=True)
            outname = os.path.join(outdir, base + "_plot.png")
        else:
            outname = os.path.splitext(filename)[0] + "_plot.png"

        plt.savefig(outname, dpi=300)
        plt.close()

    except Exception as e:
        print(f"Error processing {filename}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="MSD text files to analyze")
    parser.add_argument("--outdir", default=None,
                        help="Optional directory to save PNG plots")
    args = parser.parse_args()

    for file in args.files:
        analyze_msd(file, outdir=args.outdir)
