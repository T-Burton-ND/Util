import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.lib.distances import distance_array, calc_bonds

# ===============================
# Settings
# ===============================
data_file = "equilibrated.data"
traj_file = "production.lammpstrj"

dt_fs = 5000.0  # frame spacing in fs
rdf_range = (0.0, 12.0)
nbins = 200
li_type = 2
cl_type = 1
mass_tol = 0.2

# --- color palette ---
color_list = [
    (0.05, 0.35, 0.75), (0.9, 0.3, 0.05),
    (0.05, 0.8, 0.6), (0.9, 0.5, 0.7)
]

# ===============================
# Core Computations
# ===============================
def load_universe():
    if not os.path.exists(data_file) or not os.path.exists(traj_file):
        raise FileNotFoundError("Missing equilibrated.data or production.lammpstrj in current directory.")
    u = mda.Universe(data_file, traj_file, topology_format="DATA", format="LAMMPSDUMP", dt=dt_fs)
    return u


def compute_li_cl_rdf(u):
    li = u.select_atoms(f"type {li_type}")
    cl = u.select_atoms(f"type {cl_type}")
    rdf = InterRDF(li, cl, nbins=nbins, range=rdf_range)
    rdf.run()
    bins, g_r = rdf.results.bins, rdf.results.rdf

    box = u.dimensions[:3]
    V = np.prod(box)
    rho = len(li) / V
    cn = 4 * np.pi * rho * cumulative_trapezoid(bins**2 * g_r, bins, initial=0)

    df = pd.DataFrame({"r(Å)": bins, "g(r)": g_r, "N(r)": cn})
    df.to_csv("LiCl_RDF.csv", index=False)
    print(" Saved LiCl_RDF.csv")
    return bins, g_r, cn


def compute_li_o_rdf(u):
    type_mass_pairs = {(atom.type, atom.mass) for atom in u.atoms}
    oxy_types = [t for t, m in type_mass_pairs if np.isclose(m, 15.9, atol=mass_tol)]
    if not oxy_types:
        raise ValueError("No oxygen atom types found (mass ≈ 15.9)")
    li = u.select_atoms(f"type {li_type}")
    oxy_sel = " or ".join([f"type {t}" for t in oxy_types])
    oxy = u.select_atoms(oxy_sel)

    rdf = InterRDF(li, oxy, nbins=nbins, range=rdf_range)
    rdf.run()
    bins, g_r = rdf.results.bins, rdf.results.rdf

    df = pd.DataFrame({"r(Å)": bins, "g(r)": g_r})
    df.to_csv("LiO_RDF.csv", index=False)
    print(" Saved LiO_RDF.csv")
    return bins, g_r


def compute_li_li_rdf(u):
    li = u.select_atoms(f"type {li_type}")
    rdf = InterRDF(li, li, nbins=nbins, range=rdf_range, exclusion_block=(1, 1))
    rdf.run()
    bins, g_r = rdf.results.bins, rdf.results.rdf

    box = u.dimensions[:3]
    V = np.prod(box)
    rho = len(li) / V
    cn = 4 * np.pi * rho * cumulative_trapezoid(bins**2 * g_r, bins, initial=0)

    df = pd.DataFrame({"r(Å)": bins, "g(r)": g_r, "N(r)": cn})
    df.to_csv("LiLi_RDF.csv", index=False)
    print(" Saved LiLi_RDF.csv")
    return bins, g_r, cn


def compute_li_li_autocorr(u, r_min=2.0, r_max=4.0):
    li = u.select_atoms(f"type {li_type}")
    idx_li = li.indices
    N = len(u.atoms)
    T = len(u.trajectory)

    print(f"Computing Li–Li autocorrelation from {T} frames...")
    coords = np.empty((T, N, 3), dtype=np.float32)
    boxes = np.empty((T, 6), dtype=np.float32)
    for t, ts in enumerate(u.trajectory):
        coords[t] = u.atoms.positions
        boxes[t] = u.dimensions

    C_tau_acc = np.zeros(T)
    n_tau = np.zeros(T, dtype=np.int64)

    for t0 in range(T):
        box0 = boxes[t0]
        dmat = distance_array(coords[t0, idx_li], coords[t0], box=box0)
        i_rel, j = np.where((dmat >= r_min) & (dmat <= r_max))
        if i_rel.size == 0:
            continue
        i = idx_li[i_rel]
        mask = i != j
        i, j = i[mask], j[mask]
        if i.size == 0:
            continue

        max_tau = T - t0
        for tau in range(max_tau):
            box_tau = boxes[t0 + tau]
            rij = calc_bonds(coords[t0 + tau, i], coords[t0 + tau, j], box=box_tau)
            H_tau = (rij >= r_min) & (rij <= r_max)
            C_tau_acc[tau] += H_tau.mean()
            n_tau[tau] += 1

    valid = n_tau > 0
    C = np.zeros(T)
    C[valid] = C_tau_acc[valid] / n_tau[valid]
    time_ps = np.arange(T) * (u.trajectory.dt / 1000.0)

    df = pd.DataFrame({"time(ps)": time_ps[valid], "C(t)": C[valid]})
    df.to_csv("LiLi_autocorr.csv", index=False)
    print(" Saved LiLi_autocorr.csv")
    return time_ps[valid], C[valid]


# ===============================
# Main execution
# ===============================
print("\n=== Loading Universe ===")
u = load_universe()

print("\n--- Computing Li–Cl RDF ---")
bins1, g1, cn1 = compute_li_cl_rdf(u)

print("\n--- Computing Li–O RDF ---")
bins2, g2 = compute_li_o_rdf(u)

print("\n--- Computing Li–Li RDF ---")
bins3, g3, cn3 = compute_li_li_rdf(u)

print("\n--- Computing Li–Li Autocorrelation ---")
t, C = compute_li_li_autocorr(u)

print("\n All RDFs and autocorrelation computed and saved. Now plotting...\n")

# ===============================
# Plot all in one figure
# ===============================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

# --- Li–Cl ---
ax = axes[0]
color = color_list[0]
ax2 = ax.twinx()
ax.plot(bins1, g1, lw=3, color=color, label="Li–Cl g(r)")
ax2.plot(bins1, cn1, lw=2, ls="--", color=color, alpha=0.8)
ax.set_title("Li–Cl RDF + CN", fontsize=14, fontweight='bold')
ax.set_xlabel("r [Å]")
ax.set_ylabel("g(r)")
ax2.set_ylabel("N(r)")
ax.grid(alpha=0.3, ls='--')

# --- Li–O ---
ax = axes[1]
color = color_list[1]
ax.plot(bins2, g2, lw=3, color=color, label="Li–O g(r)")
ax.set_title("Li–O RDF", fontsize=14, fontweight='bold')
ax.set_xlabel("r [Å]")
ax.set_ylabel("g(r)")
ax.grid(alpha=0.3, ls='--')

# --- Li–Li ---
ax = axes[2]
color = color_list[2]
ax2 = ax.twinx()
ax.plot(bins3, g3, lw=3, color=color, label="Li–Li g(r)")
ax2.plot(bins3, cn3, lw=2, ls="--", color=color, alpha=0.8)
ax.set_title("Li–Li RDF + CN", fontsize=14, fontweight='bold')
ax.set_xlabel("r [Å]")
ax.set_ylabel("g(r)")
ax2.set_ylabel("N(r)")
ax.grid(alpha=0.3, ls='--')

# --- Autocorrelation ---
ax = axes[3]
color = color_list[3]
ax.plot(t, C, lw=3, color=color)
ax.set_title("Li–Li Contact Autocorrelation", fontsize=14, fontweight='bold')
ax.set_xlabel("Time [ps]")
ax.set_ylabel(r"$\langle h(0)h(t)\rangle$")
ax.grid(alpha=0.3, ls='--')

# --- Formatting ---
for ax in axes:
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='both', which='major', labelsize=10, width=1.5, length=5, direction='in')

plt.tight_layout()
plt.suptitle("Li-based System Analysis: RDFs, CNs, and Autocorrelation", fontsize=16, fontweight='bold', y=1.02)
plt.savefig("Li_all_analysis.png", dpi=300)

print("\n Combined plot saved as Li_all_analysis.png\n")