#!/usr/bin/env python3
"""
compute_diffusion_mdanalysis.py  (production-safe windowing)

What’s new (relevant to long productions with no step reset):
- Time is zeroed at the first trajectory frame by default (t=0 at production start).
- Optional absolute-step windowing: --fit-start-step/--fit-end-step
- Optional expectation check: --expect-first-step N (warns if mismatch)
- Prints a sanity banner (first/last step, total frames, ps/ns span, chosen window)
"""

from __future__ import annotations
import argparse, gzip, io, os
from pathlib import Path
from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import msd as mda_msd
import warnings
warnings.filterwarnings("ignore", message="Reader has no dt information, set to 1.0 ps")
warnings.filterwarnings("ignore", message="Guessed all Masses to 1.0")

# ----------------------------- Small Utilities ----------------------------- #

def open_maybe_gz(path: Path) -> io.TextIOBase:
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r")

def read_lammpstrj_timesteps(trj_path: Path) -> np.ndarray:
    steps: List[int] = []
    with open_maybe_gz(trj_path) as fh:
        for line in fh:
            if line.startswith("ITEM: TIMESTEP"):
                ts = fh.readline()
                if not ts:
                    break
                steps.append(int(ts.strip()))
    arr = np.asarray(steps, dtype=int)
    if arr.size > 1 and not np.all(arr[1:] >= arr[:-1]):
        raise ValueError("TIMESTEPs are not non-decreasing; check file order.")
    return arr

def compute_step_stride(steps: np.ndarray) -> Optional[int]:
    if steps.size < 3: return None
    diffs = np.diff(steps)
    # Use median to be robust to occasional missing frames
    return int(np.median(diffs))

def build_time_ps_from_steps(steps: np.ndarray, dt_fs: float, zero_time: bool) -> np.ndarray:
    t_ps_abs = steps.astype(float) * (dt_fs / 1000.0)   # fs -> ps
    return (t_ps_abs - t_ps_abs[0]) if zero_time else t_ps_abs

def pick_fit_indices(times_ps: np.ndarray,
                     tmin_ps: Optional[float],
                     tmax_ps: Optional[float]) -> Tuple[int, int]:
    if times_ps.size < 5:
        return 0, times_ps.size
    t0, t1 = times_ps[0], times_ps[-1]
    if tmin_ps is None: tmin_ps = t0 + 0.50*(t1 - t0)
    if tmax_ps is None: tmax_ps = t0 + 0.90*(t1 - t0)
    if tmax_ps <= tmin_ps:
        tmin_ps = t0 + 0.50*(t1 - t0); tmax_ps = t1
    i0 = int(np.searchsorted(times_ps, tmin_ps, side="left"))
    i1 = int(np.searchsorted(times_ps, tmax_ps, side="right"))
    i0 = max(0, min(i0, times_ps.size - 2))
    i1 = max(i0 + 2, min(i1, times_ps.size))
    return i0, i1

def polyfit_slope(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2: raise ValueError("Not enough points for a linear fit.")
    m, _b = np.polyfit(x, y, 1)
    return float(m)

def slope_to_D_cm2s(slope_A2_per_ps: float) -> float:
    return slope_A2_per_ps * (1e-4 / 6.0)

def ensure_outdir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

# ----------------------------- Core Calculator ----------------------------- #
def compute_msd_with_mdanalysis(universe: mda.Universe, selection: str,
                                stop_index: Optional[int] = None) -> np.ndarray:
    ag = universe.select_atoms(selection)
    if ag.n_atoms == 0:
        raise ValueError(f"Selection '{selection}' matched 0 atoms.")
    E = mda_msd.EinsteinMSD(universe, select=selection, msd_type="xyz", unwrap=False)
    E.run(stop=stop_index)  # <-- limit frames actually analyzed
    msd = getattr(getattr(E, "results", E), "msd", getattr(getattr(E, "results", E), "timeseries", None))
    if msd is None:
        raise AttributeError("EinsteinMSD results missing; check MDAnalysis version.")
    return np.asarray(msd, dtype=float)



def fit_diffusion(times_ps: np.ndarray, msd_A2: np.ndarray,
                  tmin_ps: Optional[float], tmax_ps: Optional[float]) -> Tuple[float, float, Tuple[int,int]]:
    n = min(times_ps.size, msd_A2.size)
    times_ps = times_ps[:n]; msd_A2 = msd_A2[:n]
    i0, i1 = pick_fit_indices(times_ps, tmin_ps, tmax_ps)
    slope = polyfit_slope(times_ps[i0:i1], msd_A2[i0:i1])
    D = slope_to_D_cm2s(slope)
    return slope, D, (i0, i1)

def plot_msd_with_fit(times_ps: np.ndarray, msd_A2: np.ndarray,
                      slope: float, i0: int, out_png: Path, title: str) -> None:
    n = min(times_ps.size, msd_A2.size)
    times_ps = times_ps[:n]; msd_A2 = msd_A2[:n]
    y0 = msd_A2[i0] - slope*times_ps[i0]
    fit_line = slope*times_ps + y0
    plt.figure()
    plt.plot(times_ps, msd_A2, lw=2, label="MSD (MDAnalysis)")
    plt.plot(times_ps, fit_line, "--", label="Linear fit")
    plt.xlabel("Time (ps)")
    plt.ylabel(r"MSD ($\mathrm{\AA^2}$)")
    plt.title(title)
    plt.legend(); plt.tight_layout()
    plt.savefig(out_png, dpi=200); plt.close()

def save_msd_csv(times_ps: np.ndarray, msd_A2: np.ndarray, out_csv: Path) -> None:
    n = min(times_ps.size, msd_A2.size)
    arr = np.column_stack([times_ps[:n], msd_A2[:n]])
    np.savetxt(out_csv, arr, delimiter=",", header="time_ps,MSD_A2", comments="")

# ----------------------------------- Main ---------------------------------- #

def run(trj_path: Path, dt_fs: float, out_dir: Path,
        tmin_ps: Optional[float], tmax_ps: Optional[float],
        zero_time: bool,
        fit_start_step: Optional[int], fit_end_step: Optional[int],
        expect_first_step: Optional[int],
        max_ps: Optional[float] = None) -> None:

    ensure_outdir(out_dir)

    # Load steps and build time axis (zeroed by default)
    steps = read_lammpstrj_timesteps(trj_path)
    if steps.size == 0:
        raise ValueError("No TIMESTEP records found in trajectory.")
    if expect_first_step is not None and steps[0] != expect_first_step:
        print(f"[warn] First frame step={steps[0]} != expected {expect_first_step}.")

    stride = compute_step_stride(steps)
    times_ps = build_time_ps_from_steps(steps, dt_fs, zero_time=zero_time)
    # --- after building times_ps ---
    i_cut = None
    if max_ps is not None:
        cutoff_ps = (0.0 if zero_time else times_ps[0]) + max_ps
        i_cut = int(np.searchsorted(times_ps, cutoff_ps, side="right"))
        if i_cut < times_ps.size:
            steps = steps[:i_cut]
            times_ps = times_ps[:i_cut]
            print(f"[note] Truncated to first {max_ps:.1f} ps → {len(times_ps)} frames retained.")

# --- (keep your early clip if you want) ---

# If user specified absolute-step window, convert to ps (respecting zero_time)
    if fit_start_step is not None or fit_end_step is not None:
        s0 = steps[0]
        def step_to_ps(S):
            if S is None: return None
            base = (S - (0 if zero_time else s0))
            return base * (dt_fs / 1000.0)
        tmin_ps = step_to_ps(fit_start_step)
        tmax_ps = step_to_ps(fit_end_step)

# >>> clip again AFTER conversion (and after truncation) <<<
    if tmin_ps is not None or tmax_ps is not None:
        t_lo, t_hi = times_ps[0], times_ps[-1]
        if tmin_ps is not None and tmin_ps < t_lo: tmin_ps = t_lo
        if tmax_ps is not None and tmax_ps > t_hi: tmax_ps = t_hi
        if tmin_ps is not None and tmax_ps is not None and tmax_ps <= tmin_ps:
            print("[note] Fit window collapsed after truncation; using default (50–90% of span).")
            tmin_ps = None
            tmax_ps = None


    # Sanity banner
    span_ps = times_ps[-1] - times_ps[0]
    print("\n=== Trajectory sanity ===")
    print(f"frames         : {steps.size}")
    print(f"first/last step: {steps[0]} / {steps[-1]}")
    if stride: print(f"median stride  : {stride} steps  (~{stride*dt_fs/1000:.3f} ps)")
    print(f"time span      : {span_ps:.3f} ps  ({span_ps/1000.0:.3f} ns)")
    if tmin_ps is not None or tmax_ps is not None:
        print(f"fit window req : tmin={tmin_ps} ps, tmax={tmax_ps} ps (times {'zeroed' if zero_time else 'absolute'})")
    print("=========================\n")

    # Load trajectory into MDAnalysis
    u = mda.Universe(str(trj_path), format="LAMMPSDUMP")

    print(f"[args] zero_time={zero_time} dt_fs={dt_fs} "
      f"fit_start_step={fit_start_step} fit_end_step={fit_end_step} "
      f"tmin_ps={tmin_ps} tmax_ps={tmax_ps}")
    
    species = [("Cl", "type 1"), ("Li", "type 2")]
    summary_lines = []
    for label, sel in species:
        print(f"[MDAnalysis] Computing MSD for {label} ({sel})...")
        msd = compute_msd_with_mdanalysis(u, selection=sel, stop_index=i_cut)


        # Align time vector with MSD length (robust to truncated MSD)
        n = min(len(msd), len(times_ps))
        tvec = times_ps[:n]; msd = msd[:n]

        slope, D, (i0, i1) = fit_diffusion(tvec, msd, tmin_ps, tmax_ps)
        save_msd_csv(tvec, msd, out_dir / f"msd_{label}.csv")
        plot_msd_with_fit(tvec, msd, slope, i0, out_dir / f"msd_fit_mdanalysis_{label}.png",
                          title=f"MDAnalysis MSD & fit: {label}")

        line = (f"{label}: fit {tvec[i0]:.2f}–{tvec[i1-1]:.2f} ps  →  "
                f"D = {D:.3e} cm^2/s  (slope={slope:.3e} Å^2/ps)")
        print("[Result]", line)
        summary_lines.append(line)

    (out_dir / "diffusion_summary.txt").write_text("\n".join(summary_lines) + "\n")
    print(f"\nSummary written to: {out_dir/'diffusion_summary.txt'}")

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute Li/Cl diffusion from a LAMMPS dump using MDAnalysis EinsteinMSD.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("trajectory", help="LAMMPS dump (supports .gz). Must include unwrapped xu yu zu.")
    p.add_argument("--dt-fs", type=float, required=True, help="Simulation timestep (fs).")
    p.add_argument("--out", type=Path, default=Path("msd_out"), help="Output directory.")
    # Fit window in time since first frame (default) or absolute time if --no-zero-time
    p.add_argument("--tmin-ps", type=float, default=None, help="Lower bound (ps) for linear fit.")
    p.add_argument("--tmax-ps", type=float, default=None, help="Upper bound (ps) for linear fit.")
    # Absolute step-based windowing (overrides tmin/tmax if provided)
    p.add_argument("--fit-start-step", type=int, default=None, help="Fit window start at this TIMESTEP.")
    p.add_argument("--fit-end-step", type=int, default=None, help="Fit window end at this TIMESTEP.")
    # Clock control & expectation check
    p.add_argument("--no-zero-time", action="store_true",
                   help="Use absolute time (do NOT zero at first frame).")
    p.add_argument("--expect-first-step", type=int, default=None,
                   help="Warn if the first TIMESTEP in the file differs from this.")
    p.add_argument("--max-ps", type=float, default=None,
               help="Analyze only the first this many picoseconds (relative to the first frame).")

    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run(
        trj_path=Path(args.trajectory),
        dt_fs=args.dt_fs,
        out_dir=Path(args.out),
        tmin_ps=args.tmin_ps,
        tmax_ps=args.tmax_ps,
        zero_time=(not args.no_zero_time),
        fit_start_step=args.fit_start_step,
        fit_end_step=args.fit_end_step,
        expect_first_step=args.expect_first_step,
        max_ps=args.max_ps,
    )
