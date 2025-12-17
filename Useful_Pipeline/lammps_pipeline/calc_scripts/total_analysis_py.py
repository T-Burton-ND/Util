import os
import sys
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.transformations as tf
from collections import Counter
from MDAnalysis.lib.distances import calc_bonds

# ===============================================================
# Element mapping for XYZ
# ===============================================================
mass_to_element = {
    1.008: "H", 6.941: "Li", 12.011: "C",
    14.007: "N", 15.999: "O", 18.998: "F",
    35.453: "Cl", 18.039: "He", 75.00: "As"
}

def build_mass_to_element_map(masses):
    masses = np.asarray(masses)
    known = np.array(list(mass_to_element.keys()))
    elements = np.empty(masses.shape, dtype=object)

    for i, m in enumerate(masses):
        idx = np.argmin(np.abs(known - m))
        elements[i] = mass_to_element[known[idx]]
    return elements


# ===============================================================
# MAIN: folder argument
# ===============================================================
if len(sys.argv) != 2:
    print("Usage: python analyze.py <simulation_folder>")
    sys.exit(1)

sim_folder = os.path.abspath(sys.argv[1])
if not os.path.isdir(sim_folder):
    print(f"Error: {sim_folder} is not a directory.")
    sys.exit(1)

sim_name = os.path.basename(sim_folder.rstrip("/"))
outdir = sim_folder + "_analysis"
os.makedirs(outdir, exist_ok=True)

print(f"Output → {outdir}")

data_file = os.path.join(sim_folder, "equilibrated.data")
traj_file = os.path.join(sim_folder, "production.lammpstrj")

if not os.path.exists(data_file) or not os.path.exists(traj_file):
    print("Error: equilibrated.data or production.lammpstrj missing inside folder.")
    sys.exit(1)

# ===============================================================
# Settings
# ===============================================================
li_type = 1
cutoff = 3.0
dt_fs = 1000.0

# ===============================================================
# Load Universe ONCE
# ===============================================================
u = mda.Universe(
    data_file, traj_file,
    topology_format="DATA",
    format="LAMMPSDUMP",
    dt=dt_fs
)

u.trajectory.add_transformations(tf.wrap(u.atoms, compound="residues"))

li = u.select_atoms(f"type {li_type}")
if len(li) != 1:
    raise ValueError("ERROR: Expected exactly ONE Li atom.")

li_idx = li.indices[0]
N = len(u.atoms)
T = len(u.trajectory)
atom_indices = np.arange(N)

coords = np.empty((T, N, 3), dtype=np.float32)
boxes  = np.empty((T, 6), dtype=np.float32)

symbols = build_mass_to_element_map(u.atoms.masses)
u.add_TopologyAttr("names", symbols)

# ===============================================================
# Precompute molecule types
# ===============================================================
resids = u.atoms.resids
types  = u.atoms.types

mol_to_types = {
    mol_id: tuple(sorted(types[resids == mol_id]))
    for mol_id in np.unique(resids)
}

unique_sigs = sorted(set(mol_to_types.values()))
typekey_dict = {sig: i for i, sig in enumerate(unique_sigs)}
molid_to_typekey = {mol: typekey_dict[sig] for mol, sig in mol_to_types.items()}

# ===============================================================
# FIRST PASS: coords, clusters, neighbors, and XYZ writing
# ===============================================================
cluster_list = []
unique_type_count_dicts = set()
neighbors_t0 = []

xyz_ext_path = os.path.join(outdir, "trajectory_ext.xyz")
f_xyz = open(xyz_ext_path, "w")

for t, ts in enumerate(u.trajectory):
    pos = u.atoms.positions
    box = u.dimensions.copy()

    coords[t] = pos
    boxes[t]  = box

    # Li–all distances
    d = np.linalg.norm(pos - pos[li_idx], axis=1)
    within = np.where((d <= cutoff) & (atom_indices != li_idx))[0]

    if within.size > 0:
        neighbors_t0.append((np.full(within.size, li_idx, dtype=int), within))
    else:
        neighbors_t0.append((None, None))

    # Cluster
    mol_ids = frozenset(resids[within])
    cluster_list.append(mol_ids)

    type_keys = [molid_to_typekey[m] for m in mol_ids]
    unique_type_count_dicts.add(tuple(sorted(Counter(type_keys).items())))

    # EXTENDED XYZ output
    lx, ly, lz, _, _, _ = box
    f_xyz.write(f"{N}\n")
    f_xyz.write(
        f'Lattice="{lx} {ly} {lz} 0 0 0" '
        f'Properties=species:S:1:pos:R:3 frame={t}\n'
    )
    for name, p in zip(u.atoms.names, pos):
        f_xyz.write(f"{name} {p[0]} {p[1]} {p[2]}\n")

f_xyz.close()

print("Saved trajectory_ext.xyz")

# ===============================================================
# CLUSTER LIFETIMES
# ===============================================================
consecutive = []
if cluster_list:
    prev = cluster_list[0]
    count = 1
    for cur in cluster_list[1:]:
        if cur == prev:
            count += 1
        else:
            consecutive.append((prev, count))
            prev = cur
            count = 1
    consecutive.append((prev, count))

lifetimes = [c for _, c in consecutive]
avg_lifetime = sum(lifetimes)/len(lifetimes) if lifetimes else 0.0

with open(os.path.join(outdir, "cluster.txt"), "w") as f:
    f.write(f"Unique molecule-type clusters: {len(unique_type_count_dicts)}\n")
    f.write(f"Average cluster lifetime (frames): {avg_lifetime:.3f}\n")

print("Saved cluster.txt")

# ===============================================================
# SECOND PASS: AUTOCORRELATION
# ===============================================================
C_tau = np.zeros(T)
n_tau = np.zeros(T, dtype=np.int64)

for t0 in range(T):
    li_i, x_j = neighbors_t0[t0]
    if li_i is None:
        continue

    max_tau = T - t0
    for tau in range(max_tau):
        rij = calc_bonds(
            coords[t0 + tau, li_i],
            coords[t0 + tau, x_j],
            box=boxes[t0 + tau]
        )
        mask = (rij <= cutoff)
        C_tau[tau] += mask.mean()
        n_tau[tau] += 1

valid = n_tau > 0
C_tau[valid] /= n_tau[valid]

time_ps = np.arange(T) * (dt_fs/1000)

df = pd.DataFrame({"time(ps)": time_ps[valid], "C(t)": C_tau[valid]})
df.to_csv(os.path.join(outdir, "Li_all_within3A_autocorr.csv"), index=False)

print("Saved Li_all_within3A_autocorr.csv")
print("\n--- ALL ANALYSIS COMPLETE ---\n")