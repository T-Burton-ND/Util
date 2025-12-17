import numpy as np
from ase.io import iread
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import argparse

# === Parameters ===
cutoff_shell = 4.0          # Å cutoff for solvation shell
TARGET_SYMBOL = "Li"        # central atom type to analyze shells around
EXCLUDE = {"Cl", "Li"}      # elements to exclude from shell analysis (e.g. {"Cl","Li"})
MIN_MOLECULES = 3           # minimum number of molecules required in a shell

# ------------------------------------------------------------
# Step 1: Read data file
# ------------------------------------------------------------
def read_data_file(datafile):
    """Read atom-to-molecule mapping and type→element mapping from LAMMPS data file."""
    atom_to_mol = {}
    mol_to_atoms = defaultdict(list)
    type_to_element = {}
    in_atoms, in_masses = False, False

    with open(datafile, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Start/stop Masses section
            if line.lower().startswith("masses"):
                in_masses = True
                continue
            if in_masses and line[0].isalpha():
                in_masses = False
            if in_masses:
                parts = line.split()
                if len(parts) >= 2:
                    type_id = int(parts[0])
                    mass = float(parts[1])
                    # crude mapping via known masses
                    mass_map = {
                        1.008: "H", 6.941: "Li", 12.011: "C",
                        14.007: "N", 15.999: "O", 18.998: "F", 35.453: "Cl"
                    }
                    closest = min(mass_map, key=lambda x: abs(x - mass))
                    type_to_element[type_id] = mass_map[closest]

            # Start/stop Atoms section
            if line.lower().startswith("atoms"):
                in_atoms = True
                continue
            if in_atoms and line[0].isalpha():
                in_atoms = False
            if in_atoms:
                parts = line.split()
                if len(parts) >= 3:
                    atom_id = int(parts[0]) - 1
                    mol_id = int(parts[1])
                    atom_type = int(parts[2])
                    atom_to_mol[atom_id] = (mol_id, atom_type)
                    mol_to_atoms[mol_id].append(atom_type)

    return atom_to_mol, mol_to_atoms, type_to_element

# ------------------------------------------------------------
# Step 2: Build molecule signatures
# ------------------------------------------------------------
def build_molecule_signatures(mol_to_atoms):
    """Return mol_id → unique signature (tuple of atom types)."""
    mol_signatures = {}
    for mol_id, atom_types in mol_to_atoms.items():
        signature = tuple(sorted(atom_types))
        mol_signatures[mol_id] = signature
    return mol_signatures

# ------------------------------------------------------------
# Step 3: Frame analyzer
# ------------------------------------------------------------
def analyze_frame(atoms, atom_to_mol, mol_signatures, type_to_element, cutoff_shell=3.0):
    """Analyze one frame and return solvation shell counter."""

    # Convert type IDs (ASE "numbers") to chemical symbols
    type_ids = atoms.get_array("numbers")
    symbols = [type_to_element[int(t)] for t in type_ids]
    atoms.set_chemical_symbols(symbols)

    positions = atoms.get_positions()
    shell_counter = Counter()

    target_indices = [i for i, s in enumerate(symbols) if s == TARGET_SYMBOL]

    for li_idx in target_indices:
        li_pos = positions[li_idx]
        dists = np.linalg.norm(positions - li_pos, axis=1)
        in_shell = np.where((dists < cutoff_shell) & (dists > 0))[0]

        shell_counts = Counter()
        for neigh_idx in in_shell:
            # Skip excluded elements
            if symbols[neigh_idx] in EXCLUDE:
                continue

            mol_id, _ = atom_to_mol[neigh_idx]
            signature = mol_signatures.get(mol_id)
            if signature:
                shell_counts[signature] += 1

        # Only keep shells with ≥ MIN_MOLECULES total
        if shell_counts and sum(shell_counts.values()) >= MIN_MOLECULES:
            frozen = frozenset(shell_counts.items())
            shell_counter[frozen] += 1

    return shell_counter

# ------------------------------------------------------------
# Step 4: Multiprocessing wrapper
# ------------------------------------------------------------
def process_frames(traj_file, atom_to_mol, mol_signatures, type_to_element, cutoff_shell=3.0):
    global_shell_counter = Counter()
    frames = list(iread(traj_file))  # load all frames into list

    args = [(frame, atom_to_mol, mol_signatures, type_to_element, cutoff_shell) for frame in frames]

    with Pool(cpu_count()) as pool:
        for result in pool.starmap(analyze_frame, args):
            global_shell_counter.update(result)

    return global_shell_counter

# ------------------------------------------------------------
# Step 5: Main script
# ------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unique solvation shell analysis from trajectory (parallel).")
    parser.add_argument("run", help="Prefix for input files (e.g., 'run_25_1')")
    args = parser.parse_args()

    run = args.run
    data_file = f"equilibrated.data"
    traj_file = f"production.lammpstrj"

    #print(f"Reading data file {data_file} ...")
    atom_to_mol, mol_to_atoms, type_to_element = read_data_file(data_file)
    mol_signatures = build_molecule_signatures(mol_to_atoms)

    # Report molecules
    all_signatures = [tuple(sorted(atom_types)) for atom_types in mol_to_atoms.values()]
    unique_signatures = set(all_signatures)
    #print(f"   → Total molecules in system: {len(mol_to_atoms)}")
    #print(f"   → Unique molecule types: {len(unique_signatures)}")

    ##print(f"Analyzing trajectory {traj_file} in parallel using {cpu_count()} cores ...")
    global_shell_counter = process_frames(traj_file, atom_to_mol, mol_signatures, type_to_element, cutoff_shell=cutoff_shell)

    #print("Analysis complete!")
    print(f"{run} {len(global_shell_counter)}")