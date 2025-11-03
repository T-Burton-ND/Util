#!/usr/bin/env python3
"""
quick_react.py — Rocking GIFs from XYZ trajectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Loads a multi-frame XYZ trajectory, infers bonds manually using covalent
radii and a bond cutoff, and creates a rocking (forward-reverse) GIF
animation using matplotlib. Works with older versions of ASE.

Usage
-----
python quick_react.py your_trajectory.xyz [OPTIONS]

Options
-------
--bond-cutoff FLOAT    Extra Ångström tolerance added to covalent radii (default: 0.45)
--size PX PX           Figure size in pixels W H (default: 400 400)
--dpi INT              Matplotlib DPI (default: 100)
--ms INT               Frame duration in milliseconds (default: 200)
--frame-digits INT     Zero padding for frame filenames (default: 3)
--outdir PATH          Folder for PNG frames (default: <stem>_frames)
--gif PATH             Output GIF (default: <stem>_reaction.gif)

Examples
--------
# Default settings
python quick_react.py reaction_path.xyz

# Larger frames and slower animation
python quick_react.py reaction_path.xyz --size 600 600 --ms 300

# Tighter bonds and custom output names/locations
python quick_react.py reaction_path.xyz --bond-cutoff 0.3 --outdir ./frames --gif out.gif

Dependencies
------------
ase, matplotlib, pillow, numpy
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from PIL import Image

import matplotlib.pyplot as plt
from ase.io import read
from ase.geometry import get_distances
from ase.data import covalent_radii
from ase.data.colors import cpk_colors

__version__ = "0.2.0"

# ──────────────────────────────────────────────────────────────
#                   Helper: Draw Atoms with Bonds
# ──────────────────────────────────────────────────────────────
def draw_atoms_with_bonds(atoms, ax, *, bond_cutoff: float = 0.45, padding: float = 1.5) -> None:
    """
    Draw atoms and bonds in 2D (x-y projection).
    Automatically adjusts view to prevent cutoff at edges.
    """
    positions = atoms.get_positions()
    numbers = atoms.get_atomic_numbers()
    xy = positions[:, :2]
    radii = covalent_radii[numbers] * 0.4
    colors = [cpk_colors[Z] if Z < len(cpk_colors) else (0.5, 0.5, 0.5) for Z in numbers]

    # Draw atoms
    for (x, y), r, c in zip(xy, radii, colors):
        circle = plt.Circle((x, y), r, color=c, ec="black", lw=0.3, zorder=2)
        ax.add_patch(circle)

    # Draw bonds
    cutoff_matrix = (
        covalent_radii[numbers][:, None] + covalent_radii[numbers][None, :] + bond_cutoff
    )
    raw_distances, _ = get_distances(positions, positions)
    distances = np.linalg.norm(raw_distances, axis=-1) if raw_distances.shape[-1] == 3 else raw_distances

    bonded = (distances < cutoff_matrix) & (distances > 0.1)
    n = len(atoms)
    for i in range(n):
        for j in range(i):
            if bonded[i, j]:
                x1, y1 = xy[i]
                x2, y2 = xy[j]
                ax.plot([x1, x2], [y1, y2], color="black", linewidth=1.0, zorder=1)

    # Adjust limits to fit everything with padding
    x_vals, y_vals = xy[:, 0], xy[:, 1]
    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()
    x_pad = (x_max - x_min) * padding * 0.5 if x_max > x_min else 1.0
    y_pad = (y_max - y_min) * padding * 0.5 if y_max > y_min else 1.0
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_aspect("equal")
    ax.axis("off")


# ──────────────────────────────────────────────────────────────
#                         CLI / Main
# ──────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="quick_react.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Create a rocking (forward-then-back) GIF from a multi-frame XYZ trajectory.",
        epilog=__doc__.split("Examples", 1)[1],  # reuse Examples block from top docstring
    )
    p.add_argument("xyz", help="Path to multi-frame XYZ trajectory")
    p.add_argument("--bond-cutoff", type=float, default=0.45, help="Extra Å tolerance added to covalent radii")
    p.add_argument("--size", nargs=2, type=int, metavar=("W", "H"), default=(400, 400), help="Figure size in pixels")
    p.add_argument("--dpi", type=int, default=100, help="Matplotlib DPI")
    p.add_argument("--ms", type=int, default=200, help="Frame duration in milliseconds")
    p.add_argument("--frame-digits", type=int, default=3, help="Zero padding for frame filenames")
    p.add_argument("--outdir", type=str, help="Folder for PNG frames (default: <stem>_frames)")
    p.add_argument("--gif", type=str, help="Output GIF path (default: <stem>_reaction.gif)")
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    input_xyz = args.xyz
    if not os.path.isfile(input_xyz):
        print(f"Error: file not found: {input_xyz}", file=sys.stderr)
        return 2

    stem = Path(input_xyz).stem
    outdir = args.outdir or f"{stem}_frames"
    gif_output = args.gif or f"{stem}_reaction.gif"
    image_prefix = "frame"
    frame_digits = args.frame_digits
    dpi = args.dpi
    image_size = tuple(args.size)
    bond_cutoff = float(args.bond_cutoff)
    frame_duration_ms = int(args.ms)

    # --- Load XYZ file ---
    print(f"[✓] Loading: {input_xyz}")
    try:
        atoms_list = read(input_xyz, index=":")
    except Exception as exc:
        print(f"[✗] Failed to read XYZ: {exc}", file=sys.stderr)
        return 1
    if not atoms_list:
        print("[✗] No frames found in trajectory.", file=sys.stderr)
        return 1

    # --- Rocking: forward + backward without duplicate ---
    frames = atoms_list + atoms_list[-2::-1]

    # --- Create output dir ---
    os.makedirs(outdir, exist_ok=True)

    # --- Render each frame ---
    print(f"[✓] Rendering {len(frames)} frames → {outdir}")
    for i, atoms in enumerate(frames):
        fig, ax = plt.subplots(figsize=(image_size[0] / dpi, image_size[1] / dpi), dpi=dpi)
        ax.axis("off")
        draw_atoms_with_bonds(atoms, ax, bond_cutoff=bond_cutoff)

        frame_path = os.path.join(outdir, f"{image_prefix}_{i:0{frame_digits}}.png")
        plt.savefig(frame_path, bbox_inches="tight", pad_inches=0)
        plt.close(fig)

    # --- Assemble into GIF ---
    print(f"[✓] Assembling GIF: {gif_output}")
    pngs = [f for f in sorted(os.listdir(outdir)) if f.endswith(".png")]
    if not pngs:
        print("[✗] No PNG frames found to assemble.", file=sys.stderr)
        return 1

    images = [Image.open(os.path.join(outdir, f)) for f in pngs]
    images[0].save(
        gif_output,
        save_all=True,
        append_images=images[1:],
        duration=frame_duration_ms,
        loop=0,
    )

    print(f"[✓] Done! GIF saved to: {gif_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
