#!/usr/bin/env python3
"""
rdkit_molecule_viz.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author  : Thomas J. Burton – Savoie Research Group, UND
Updated : 2025-06-27
License : MIT

Description
-----------
Generate a 2D depiction of a molecule from a SMILES string using RDKit,
with ACS 1996 element coloring by default.

Usage
-----
# Default usage (ACS 1996 colors, output: CCO_viz.png)
python rdkit_molecule_viz.py "CCO"

# Override to black & white rendering
python rdkit_molecule_viz.py "CCO" --no-multicolor

# Custom output filename
python rdkit_molecule_viz.py "CCO" -o ethanol.png

Dependencies
------------
conda install -c conda-forge rdkit
"""

import argparse
import re
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DCairo

# ------------------------------------------------------------------------------
# ACS 1996 element color palette (RGB values normalized to [0, 1])
# ------------------------------------------------------------------------------
ACS_1996_COLORS = {
    'H':  (1.0, 1.0, 1.0),     # White
    'C':  (0.0, 0.0, 0.0),     # Black
    'N':  (0.0, 0.0, 1.0),     # Blue
    'O':  (1.0, 0.0, 0.0),     # Red
    'F':  (0.0, 1.0, 0.0),     # Green
    'Cl': (0.0, 1.0, 0.0),     # Green
    'Br': (0.65, 0.16, 0.16),  # Brown
    'I':  (0.58, 0.0, 0.82),   # Purple
    'S':  (1.0, 1.0, 0.0),     # Yellow
    'P':  (1.0, 0.5, 0.0),     # Orange
}

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

def sanitize_filename(text: str) -> str:
    """Sanitize SMILES to use as a safe filename."""
    return re.sub(r'[^a-zA-Z0-9]', '_', text)

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Draw a 2D depiction of a molecule from a SMILES string."
    )
    parser.add_argument("smiles", help="Input SMILES string for the molecule")
    parser.add_argument(
        "-o", "--output",
        help="Output image file name (default: <smiles>_viz.png)"
    )
    parser.add_argument(
        "--no-multicolor", action="store_true",
        help="Disable ACS 1996 elemental coloring (use black & white)"
    )
    return parser.parse_args()

# ------------------------------------------------------------------------------
# Main drawing logic
# ------------------------------------------------------------------------------

def draw_molecule(smiles: str, output_file: str, multicolor: bool = True):
    """
    Render and save a 2D image of a molecule from SMILES input.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    output_file : str
        Path to save the .png file.
    multicolor : bool
        Use ACS-style elemental coloring if True; black and white if False.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    # Compute 2D coordinates for proper layout
    rdDepictor.Compute2DCoords(mol)

    # Set up the RDKit drawing canvas
    drawer = MolDraw2DCairo(300, 300)
    opts = drawer.drawOptions()
    opts.setBackgroundColour((1, 1, 1, 0))  # RGBA: white with 0 alpha (fully transparent)


    # Apply ACS color scheme or black & white
    if multicolor:
        palette = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in ACS_1996_COLORS:
                palette[atom.GetIdx()] = ACS_1996_COLORS[symbol]
        opts.atomColours = palette
    else:
        opts.useBWAtomPalette()

    # Draw and save
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    with open(output_file, "wb") as f:
        f.write(drawer.GetDrawingText())

    print(f"[✓] Molecule image saved as: {output_file}")

# ------------------------------------------------------------------------------
# Entrypoint
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    smiles = args.smiles.strip()

    # Determine output filename
    if args.output:
        output_path = args.output
    else:
        base = sanitize_filename(smiles)
        output_path = f"{base}_viz.png"

    draw_molecule(smiles, output_path, multicolor=not args.no_multicolor)
