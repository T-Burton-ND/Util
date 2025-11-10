#!/usr/bin/env python3
# extract_enthalpy_simple.py
import re, sys, collections

HARTREE_TO_KCAL = 627.509474

# --- fill what you have; H is a reasonable 298 K estimate (≈ -0.5 + RT) ---
ATOM_H_EH = {"H": -0.49905579, "C": -37.79035904, "O": -74.96898166}

def find_last(pattern, text):
    m = None
    for m in re.finditer(pattern, text, re.S):
        pass
    return m

def get_total_enthalpy_Eh(txt):
    m = find_last(r"Total Enthalpy\s+.*?(-?\d+\.\d+(?:[Ee][+-]?\d+)?)\s*Eh", txt)
    if m: return float(m.group(1))
    # fallback: Total Energy + Thermal Enthalpy correction
    e = find_last(r"Total Energy\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)\s*Eh", txt)
    hc = find_last(r"Thermal Enthalpy correction\s+.*?(-?\d+\.\d+(?:[Ee][+-]?\d+)?)\s*Eh", txt)
    if e and hc: return float(e.group(1)) + float(hc.group(1))
    sys.exit("Could not find enthalpy terms. Did you run with 'Freq'?")

def get_element_counts(txt):
    # Scan the last CARTESIAN COORDINATES block
    blk = find_last(r"CARTESIAN COORDINATES\s*\(ANGSTROEM\)\s*-+\s*(.*?)\n\s*-+", txt)
    counts = collections.Counter()
    if blk:
        for line in blk.group(1).splitlines():
            p = line.split()
            if len(p) >= 4 and p[0].isalpha():
                counts[p[0]] += 1
    return counts

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "job.out"
    with open(path, "r", errors="ignore") as f:
        txt = f.read()

    H_Eh = get_total_enthalpy_Eh(txt)
    H_kcal = H_Eh * HARTREE_TO_KCAL
    #print(f"Total Enthalpy: {H_Eh:.8f} Eh  ({H_kcal:.3f} kcal/mol)")

    counts = get_element_counts(txt)
    if counts:
        corr_Eh = sum(ATOM_H_EH.get(el, 0.0) * n for el, n in counts.items())
        corr_kcal = corr_Eh * HARTREE_TO_KCAL
        #print("Composition:", " ".join(f"{el}:{n}" for el, n in sorted(counts.items())))
        #print(f"Sum atomic enthalpies: {corr_Eh:.8f} Eh  ({corr_kcal:.3f} kcal/mol)")
        print(f"{H_kcal - corr_kcal:.3f}")
