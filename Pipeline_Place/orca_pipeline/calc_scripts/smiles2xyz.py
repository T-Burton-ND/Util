#!/usr/bin/env python3
import argparse, sys, shutil, subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def embed_rdkit(mol):
    ps = AllChem.ETKDGv3(); ps.useSmallRingTorsions = True; ps.randomSeed = 0xBEEF
    if AllChem.EmbedMolecule(mol, ps) == 0: return True
    mol.RemoveAllConformers()
    ps.useRandomCoords = True; ps.maxAttempts = 1000
    if AllChem.EmbedMolecule(mol, ps) == 0: return True
    mol.RemoveAllConformers()
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=ps, pruneRmsThresh=0.05)
    return len(ids) > 0

def smiles2xyz(smiles, outxyz):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Bad SMILES")
    mol = Chem.AddHs(mol)
    if not embed_rdkit(mol): return False
    try:
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s", maxIters=300)
        best = min(range(len(res)), key=lambda i: res[i][1])
        conf = mol.GetConformer(best)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500); conf = mol.GetConformer()
    with open(outxyz, "w") as f:
        n = mol.GetNumAtoms(); f.write(f"{n}\nSMILES: {smiles}\n")
        for a in mol.GetAtoms():
            p = conf.GetAtomPosition(a.GetIdx())
            f.write(f"{a.GetSymbol():2s} {p.x:.6f} {p.y:.6f} {p.z:.6f}\n")
    return True

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles"); ap.add_argument("-o","--out", default="molecule.xyz")
    args = ap.parse_args()
    try:
        ok = smiles2xyz(args.smiles, args.out)
        if not ok and shutil.which("obabel"):
            # last-chance: Open Babel
            cmd = ["obabel", f"-:{args.smiles}", "--gen3d", "-O", args.out]
            subprocess.run(cmd, check=True)
            ok = True
        if not ok: sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"[smiles2xyz ERROR] {e}\n"); sys.exit(1)
