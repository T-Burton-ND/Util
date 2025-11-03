import numpy as np
import MDAnalysis as mda
from collections import defaultdict

# ----------------------------
# User parameters
# ----------------------------
run = "1_1"       # base name of files
r_cut = 3.0         # cutoff radius (Å)
li_type = 2         # atom type for Li
dt = 5000           # timestep between frames
smooth_window = 3   # require change to persist ≥ this many frames before accepting it
# ----------------------------

# --- Load trajectory ---
print(f"Loading trajectory: equilibrated.data / production.lammpstrj")
u = mda.Universe(f"equilibrated.data", f"production.lammpstrj",
                 topology_format="DATA", format="LAMMPSDUMP", dt=dt)

# --- Define atom selections ---
li_atoms = u.select_atoms(f"type {li_type}")
neighbors = u.select_atoms(f"not type {li_type}")
print(f"Found {len(li_atoms)} Li atoms and {len(neighbors)} neighbor atoms")

# --- Storage ---
shell_persistence = defaultdict(list)  # li_id -> list of (set_of_neighbor_ids, duration_in_frames)

# --- Iterate over Li atoms individually ---
print("Analyzing solvation shell persistence with smoothing...")
for li in li_atoms:
    prev_shell = None
    count = 0
    pending_shell = None
    pending_count = 0

    for ts in u.trajectory:
        # compute distances from this Li to all neighbor atoms
        dists = np.linalg.norm(neighbors.positions - li.position, axis=1)
        current_shell = frozenset(neighbors.ids[dists < r_cut])

        # if same as current stable shell, continue
        if current_shell == prev_shell:
            count += 1
            # reset pending change
            pending_shell = None
            pending_count = 0

        else:
            # shell appears to have changed
            if pending_shell is None:
                # start counting a possible new shell
                pending_shell = current_shell
                pending_count = 1
            elif current_shell == pending_shell:
                pending_count += 1
            else:
                # changed again before stabilization → reset pending
                pending_shell = current_shell
                pending_count = 1

            # only commit to change if new shell persisted long enough
            if pending_count >= smooth_window:
                if prev_shell is not None:
                    shell_persistence[li.id].append((prev_shell, count))
                prev_shell = pending_shell
                count = pending_count
                pending_shell = None
                pending_count = 0

    # store final shell
    if prev_shell is not None:
        shell_persistence[li.id].append((prev_shell, count))

# --- Global statistics ---
all_lifetimes = []
for li, shells in shell_persistence.items():
    lifetimes = [dur for _, dur in shells]
    all_lifetimes.extend(lifetimes)

print("\n---------------- Results ----------------")
if all_lifetimes:
    mean_lifetime = np.mean(all_lifetimes)
    std_lifetime = np.std(all_lifetimes)
    total_shells = len(all_lifetimes)
    print(f"Global mean solvation-shell lifetime : {mean_lifetime:.2f} frames ± {std_lifetime:.2f}")
    #print(f"Total number of shell events          : {total_shells}")
    #print(f"Average per Li                        : {total_shells / len(li_atoms):.1f} events/Li\n")

    # Optional: print per-Li summary
    #for li, shells in shell_persistence.items():
    #   lifetimes = [dur for _, dur in shells]
    #    print(f"Li {li}: mean lifetime = {np.mean(lifetimes):.1f} frames, "
    #          f"unique shells = {len(shells)}")
else:
    print("No solvation shell data found.")
