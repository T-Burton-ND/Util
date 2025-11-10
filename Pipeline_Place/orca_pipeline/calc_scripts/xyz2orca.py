#!/usr/bin/env python3
# xyz2orca.py
import argparse, os

def read_xyz(path):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    try:
        n = int(lines[0])
        coord_lines = lines[2:2+n]
    except Exception:
        coord_lines = lines[2:]  # fall back if header count is off
    atoms = []
    for l in coord_lines:
        parts = l.split()
        if len(parts) >= 4:
            atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms

def write_orca(atoms, out_path, method, basis, charge, mult, title):
    with open(out_path, "w") as f:
        f.write(f"! {method} {basis} TightSCF Opt Freq\n\n")
        f.write(f"%maxcore 2000\n\n")  # simple, safe default; adjust as needed
        f.write(f"# {title}\n\n")
        f.write(f"* xyz {charge} {mult}\n")
        for el,x,y,z in atoms:
            f.write(f"{el:<2} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("*\n")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="XYZ → ORCA input")
    p.add_argument("xyz", help="Input .xyz file")
    p.add_argument("-o","--out", default="orca.inp", help="Output .inp (default: <xyz basename>.inp)")
    p.add_argument("--method", default="B3LYP", help="DFT method (default: B3LYP)")
    p.add_argument("--basis", default="def2-SVP", help="Basis set (default: def2-SVP)")
    p.add_argument("-c","--charge", type=int, default=0, help="Total charge (default: 0)")
    p.add_argument("-m","--mult", type=int, default=1, help="Multiplicity (default: 1)")
    args = p.parse_args()

    atoms = read_xyz(args.xyz)
    if not atoms:
        raise SystemExit("No atoms parsed from XYZ.")
    out = args.out or os.path.splitext(os.path.basename(args.xyz))[0] + ".inp"
    title = f"Auto-generated from {os.path.basename(args.xyz)}"
    write_orca(atoms, out, args.method, args.basis, args.charge, args.mult, title)
    print(f"Wrote {out} ({len(atoms)} atoms)")
