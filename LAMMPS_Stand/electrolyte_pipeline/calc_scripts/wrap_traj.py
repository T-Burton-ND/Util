#!/usr/bin/env python3
"""
Wrap an unwrapped LAMMPS trajectory (xu yu zu) into the primary box using MDAnalysis.
Inputs (assumed in CWD): system.data, production.lammpstrj (or production.lammpstrj.gz)
Output: production_wrapped_<timestamp>.lammpstrj
"""

import os
from datetime import datetime
import MDAnalysis as mda
from MDAnalysis.coordinates.LAMMPSDUMP import DumpWriter

data_file = "equilibrated.data"
traj_candidates = ["production.lammpstrj", "production.lammpstrj.gz"]

if not os.path.exists(data_file):
    raise FileNotFoundError("Missing system.data in this directory.")

traj_file = next((f for f in traj_candidates if os.path.exists(f)), None)
if traj_file is None:
    raise FileNotFoundError("Missing production.lammpstrj (or .gz) in this directory.")

print("Loading trajectory...")
# Force coordinate reader to LAMMPSDUMP (fixes the error you saw)
u = mda.Universe(data_file, traj_file, format="LAMMPSDUMP")

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
out_file = f"production_wrapped.lammpstrj"

print(f"Wrapping and writing to {out_file} ...")
# Use the LAMMPS dump writer so we keep a LAMMPS-style dump
with DumpWriter(out_file, n_atoms=u.atoms.n_atoms) as W:
    # Keep molecules whole if residues/molecules are defined; else fall back to atom-wise
    for ts in u.trajectory:
        try:
            u.atoms.wrap(compound="residues")  # best if topology has molecules
        except Exception:
            u.atoms.wrap(compound="atoms")     # safe fallback

        # Write id, type, x,y,z; include 'mol' if available
        cols = ["id", "type", "x", "y", "z"]
        if getattr(u.atoms, "resnums", None) is not None or "mol" in u.atoms.atoms.names:
            # DumpWriter recognizes 'mol' when available via topology/residues
            # (MDAnalysis maps residues -> molecule IDs for LAMMPS dumps)
            cols = ["id", "mol", "type", "x", "y", "z"]
        W.write(u.atoms, columns=cols)

print(f"✅ Done! Wrote: {out_file}")
