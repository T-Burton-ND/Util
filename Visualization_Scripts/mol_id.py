#!/usr/bin/env python3
"""
mol_id.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author  : Thomas J. Burton – Savoie Research Group, UND
Updated : 2025-07-22
License : MIT

Description
-----------
Generate a 2D depiction of a molecule from a SMILES string or InChIKey using RDKit,
with optional ACS 1996 element coloring and PubChem name resolution.

Usage
-----
# Visualize a molecule from SMILES
python mol_id.py "CCO"

# Visualize using an InChIKey
python mol_id.pyy "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

# Black & white rendering
python mol_id.py "CCO" --no-multicolor

# Custom output filename
python mol_id.py "CCO" -o ethanol.png

Dependencies
------------
conda install -c conda-forge rdkit
pip install requests
"""

import argparse
import re
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
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
    """
    Convert a string into a filesystem-safe filename.
    
    Parameters
    ----------
    text : str
        Input string to sanitize (e.g., a SMILES string).
    
    Returns
    -------
    str
        Sanitized version safe for filenames.
    """
    return re.sub(r'[^a-zA-Z0-9]', '_', text)

def identify_input_type(query: str) -> str:
    """
    Identify whether input is a SMILES string or InChIKey.
    
    Parameters
    ----------
    query : str
        Input chemical identifier string.
    
    Returns
    -------
    str
        'smiles' or 'inchikey'
    
    Raises
    ------
    ValueError
        If input doesn't match known formats.
    """
    if re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", query):
        return "inchikey"
    elif Chem.MolFromSmiles(query):
        return "smiles"
    else:
        raise ValueError(f"Unrecognized input format: {query}")

def inchikey_to_smiles(inchikey: str) -> str:
    """
    Use PubChem API to resolve InChIKey to a canonical SMILES.
    
    Parameters
    ----------
    inchikey : str
        Valid InChIKey string.
    
    Returns
    -------
    str
        Canonical SMILES string.
    
    Raises
    ------
    ValueError
        If lookup fails.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON"
    response = requests.get(url, timeout=10)
    if response.status_code != 200:
        raise ValueError(f"Failed to resolve InChIKey: {inchikey}")
    return response.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]

def get_common_name_from_inchikey(inchikey: str) -> str | None:
    """
    Retrieve a common or trivial name from PubChem using InChIKey.
    
    Parameters
    ----------
    inchikey : str
        InChIKey for the compound.
    
    Returns
    -------
    str or None
        Best-matching common name, or None if not found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/synonyms/JSON"
    response = requests.get(url, timeout=10)
    if response.status_code != 200:
        return None
    try:
        synonyms = response.json()["InformationList"]["Information"][0]["Synonym"]
        for name in synonyms:
            if not any(x in name for x in ["CID", "CHEBI", "InChI", "C[C@H]"]):
                return name
        return synonyms[0]  # Fallback
    except Exception:
        return None

def get_molecule_info(query: str, verbose: bool = True) -> tuple[str, Chem.Mol, str | None]:
    """
    Convert an input string to RDKit Mol object and canonical SMILES.
    
    Parameters
    ----------
    query : str
        SMILES or InChIKey input string.
    verbose : bool
        If True, print resolution information.
    
    Returns
    -------
    tuple[str, Mol, str or None]
        Canonical SMILES, RDKit Mol, common name (if resolved).
    
    Raises
    ------
    ValueError
        If conversion or resolution fails.
    """
    query_type = identify_input_type(query)
    if query_type == "inchikey":
        smiles = inchikey_to_smiles(query)
        mol = Chem.MolFromSmiles(smiles)
        name = get_common_name_from_inchikey(query)
        if verbose:
            print(f"[✓] Resolved InChIKey to SMILES: {smiles}")
            if name:
                print(f"[✓] Common name: {name}")
    else:
        mol = Chem.MolFromSmiles(query)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {query}")
        smiles = Chem.MolToSmiles(mol, canonical=True)
        name = None
        if verbose:
            print(f"[✓] Canonical SMILES: {smiles}")
    return smiles, mol, name

def draw_molecule(smiles: str, output_file: str, multicolor: bool = True):
    """
    Render and save a 2D molecule image from a SMILES string.
    
    Parameters
    ----------
    smiles : str
        Canonical SMILES string.
    output_file : str
        Output .png filename.
    multicolor : bool
        Use ACS-style coloring if True; else black-and-white.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    rdDepictor.Compute2DCoords(mol)

    drawer = MolDraw2DCairo(300, 300)
    opts = drawer.drawOptions()
    opts.setBackgroundColour((1, 1, 1, 0))  # Transparent background

    if multicolor:
        palette = {
            atom.GetIdx(): ACS_1996_COLORS.get(atom.GetSymbol(), (0.0, 0.0, 0.0))
            for atom in mol.GetAtoms()
        }
        opts.atomColours = palette
    else:
        opts.useBWAtomPalette()

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    with open(output_file, "wb") as f:
        f.write(drawer.GetDrawingText())
    
    print(f"[✓] Molecule image saved as: {output_file}")

def parse_args():
    """
    Parse command-line arguments for script usage.
    
    Returns
    -------
    argparse.Namespace
        Parsed arguments including input string and output settings.
    """
    parser = argparse.ArgumentParser(
        description="Draw a 2D depiction of a molecule from a SMILES string or InChIKey."
    )
    parser.add_argument("query", help="SMILES string or InChIKey")
    parser.add_argument(
        "-o", "--output", help="Output image file name (default: <smiles>_viz.png)"
    )
    parser.add_argument(
        "--no-multicolor", action="store_true",
        help="Disable ACS 1996 elemental coloring (use black & white)"
    )
    return parser.parse_args()

# ------------------------------------------------------------------------------
# Entrypoint
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    try:
        smiles, mol, name = get_molecule_info(args.query)

        # Determine output file name
        output_file = args.output or f"{sanitize_filename(smiles)}_viz.png"

        # Draw molecule
        draw_molecule(smiles, output_file, multicolor=not args.no_multicolor)

        # Print optional name
        if name:
            print(f"[✓] Common/trivial name: {name}")

    except Exception as e:
        print(f"[✗] Error: {e}")
