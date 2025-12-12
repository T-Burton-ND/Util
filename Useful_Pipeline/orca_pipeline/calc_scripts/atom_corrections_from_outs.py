#!/usr/bin/env python3
# atom_corrections_from_outs.py
import sys, re, glob

HARTREE_TO_KCAL = 627.509474

def total_enthalpy(out_text):
    m = re.search(r"Total Enthalpy\s+(-?\d+\.\d+)\s*Eh", out_text)
    if m: return float(m.group(1))
    e = re.search(r"Total Energy\s+(-?\d+\.\d+)\s*Eh", out_text)
    hc = re.search(r"Thermal correction to Enthalpy\s+(\d+\.\d+)\s*Eh", out_text)
    if e and hc: return float(e.group(1)) + float(hc.group(1))
    raise SystemExit("Missing enthalpy terms in an atom .out")

pairs = []
for fn in sorted(glob.glob("*.out")):
    el = fn.split(".")[0]  # expects H.out, C.out, ...
    with open(fn,"r",errors="ignore") as f:
        H_Eh = total_enthalpy(f.read())
    pairs.append((el, H_Eh * HARTREE_TO_KCAL))

print("ATOM_H_CORR = {")
for el, kcal in pairs:
    print(f'    "{el}": {kcal:.6f},')
print("}")
