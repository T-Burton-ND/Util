import MDAnalysis as mda
from MDAnalysis.analysis import msd
#import matplotlib.pyplot as plt
import numpy as np
#from sklearn.linear_model import LinearRegression
#from scipy.optimize import curve_fit
import argparse

# ----------------------------                                                                                                                                                                                                                                   
# Parse command-line arguments                                                                                                                                                                                                                                   
# ----------------------------                                                                                                                                                                                                                                   
parser = argparse.ArgumentParser(description="MSD analysis with log-log least-squares fit")
parser.add_argument("run", type=str, help="Prefix for input files (expects {run}.data.lmp and {run}.diffusion.lammpstrj)")
args = parser.parse_args()
run = args.run

print(f"[MSD] Processing run: {run}")
print("[MSD] Building Universe...")
# Universe                                                                                                                                                                                                                                                       
u = mda.Universe(
    f"equilibrated.data", f"production.lammpstrj",
	topology_format="DATA", format="LAMMPSDUMP", dt=5000
)

n_atoms = u.atoms.n_atoms
n_frames = len(u.trajectory)

print(f"[MSD] Topology: equilibrated.data with {n_atoms} atoms")
print(f"[MSD] Trajectory: production.lammpstrj with {n_frames} frames")
print(f"[MSD] Frame spacing = {u.trajectory.dt:.4f} fs")

# Select only atom type 1                                                                                                                                                                                                                                        
atoms_type1 = u.select_atoms("type 2")
print(f"[MSD] Selected {atoms_type1.n_atoms} atoms of type 1")

# ----------------------------                                                                                                                                                                                                                                   
# Compute MSD with numba backend                                                                                                                                                                                                                                 
# ----------------------------                                                                                                                                                                                                                                   
print("[MSD] Computing MSD...")
MSD = msd.EinsteinMSD(atoms_type1, msd_type='xyz', fft=True, backend="numba")
MSD.run()

# Extract MSD results                                                                                                                                                                                                                                            
if "msd" in MSD.results:
    msd_vals = MSD.results["msd"]
elif "timeseries" in MSD.results:
    msd_vals = MSD.results["timeseries"]
else:
    raise KeyError(f"Could not find MSD array, keys = {MSD.results.keys()}")

# ----------------------------                                                                                                                                                                                                                                   
# Build time axis                                                                                                                                                                                                                                                
# ----------------------------                                                                                                                                                                                                                                   
nframes = len(msd_vals)
time_fs = np.arange(nframes) * u.trajectory.dt   # fs                                                                                                                                                                                                            

# ----------------------------                                                                                                                                                                                                                                   
# Write MSD to text file                                                                                                                                                                                                                                         
# ----------------------------                                                                                                                                                                                                                                   
output_file = f"mda_msd_{run}.txt"
with open(output_file, "w") as f:
    f.write("# Time_fs\tMSD_A2\n")
    for t, msd_val in zip(time_fs, msd_vals):
        f.write(f"{t:.6e}\t{msd_val:.6e}\n")
print(f"[MSD] Saved MSD data to {output_file}")

