#!/usr/bin/env python3
"""
quick_react.py

Usage:
    python quick_react.py your_trajectory.xyz

Description:
    Loads a multi-frame XYZ trajectory, infers bonds manually using covalent
    radii and a bond cutoff, and creates a rocking (forward-reverse) GIF
    animation using matplotlib. Works with older versions of ASE.
"""

import os
import sys
from pathlib import Path
from ase.io import read
from ase.visualize.plot import plot_atoms
from ase.geometry import get_distances
from ase.data import covalent_radii
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np

# ──────────────────────────────────────────────────────────────
#                   Helper: Draw Atoms with Bonds
# ──────────────────────────────────────────────────────────────
def draw_atoms_with_bonds(atoms, ax, bond_cutoff=0.45, padding=1.5):
    """
    Draw atoms and bonds in 2D (x-y projection).
    Automatically adjusts view to prevent cutoff at edges.
    """
    from ase.data import covalent_radii
    from ase.data.colors import cpk_colors

    positions = atoms.get_positions()
    numbers = atoms.get_atomic_numbers()
    xy = positions[:, :2]
    radii = covalent_radii[numbers] * 0.4
    colors = [cpk_colors[Z] if Z < len(cpk_colors) else (0.5, 0.5, 0.5) for Z in numbers]

    # Draw atoms
    for (x, y), r, c in zip(xy, radii, colors):
        circle = plt.Circle((x, y), r, color=c, ec='black', lw=0.3, zorder=2)
        ax.add_patch(circle)

    # Draw bonds
    cutoff_matrix = (
        covalent_radii[numbers][:, None] + covalent_radii[numbers][None, :] + bond_cutoff
    )
    raw_distances, _ = get_distances(positions, positions)
    if raw_distances.shape[-1] == 3:
        distances = np.linalg.norm(raw_distances, axis=-1)
    else:
        distances = raw_distances

    bonded = (distances < cutoff_matrix) & (distances > 0.1)
    for i in range(len(atoms)):
        for j in range(i):
            if bonded[i, j]:
                x1, y1 = xy[i]
                x2, y2 = xy[j]
                ax.plot([x1, x2], [y1, y2], color='black', linewidth=1.0, zorder=1)

    # Adjust limits to fit everything with padding
    x_vals, y_vals = xy[:, 0], xy[:, 1]
    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()
    x_pad = (x_max - x_min) * padding * 0.5
    y_pad = (y_max - y_min) * padding * 0.5
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_aspect('equal')
    ax.axis('off')


# ──────────────────────────────────────────────────────────────
#                         Main Script
# ──────────────────────────────────────────────────────────────

# --- Command-line input ---
if len(sys.argv) < 2:
    print("Usage: python quick_react.py your_trajectory.xyz")
    sys.exit(1)

input_xyz = sys.argv[1]
stem = Path(input_xyz).stem

# --- Parameters ---
output_dir = f"{stem}_frames"
gif_output = f"{stem}_reaction.gif"
image_prefix = "frame"
frame_digits = 3
dpi = 100
image_size = (400, 400)  # pixels
bond_cutoff = 0.45       # tolerance added to covalent radii
frame_duration_ms = 200  # duration of each frame in milliseconds
# ------------------

# --- Load XYZ file ---
print(f"[✓] Loading: {input_xyz}")
atoms_list = read(input_xyz, index=":")

# --- Rocking: forward + backward without duplicate ---
frames = atoms_list + atoms_list[-2::-1]

# --- Create output dir ---
os.makedirs(output_dir, exist_ok=True)

# --- Render each frame ---
print(f"[✓] Rendering {len(frames)} frames...")
for i, atoms in enumerate(frames):
    fig, ax = plt.subplots(figsize=(image_size[0]/dpi, image_size[1]/dpi), dpi=dpi)
    ax.axis("off")
    draw_atoms_with_bonds(atoms, ax, bond_cutoff=bond_cutoff)

    frame_path = os.path.join(output_dir, f"{image_prefix}_{i:0{frame_digits}}.png")
    plt.savefig(frame_path, bbox_inches="tight", pad_inches=0)
    plt.close()

# --- Assemble into GIF ---
print(f"[✓] Assembling GIF: {gif_output}")
images = [
    Image.open(os.path.join(output_dir, f))
    for f in sorted(os.listdir(output_dir)) if f.endswith(".png")
]

images[0].save(
    gif_output,
    save_all=True,
    append_images=images[1:],
    duration=frame_duration_ms,
    loop=0
)

print(f"[✓] Done! GIF saved to: {gif_output}")
