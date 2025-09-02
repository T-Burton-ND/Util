#!/usr/bin/env python3
"""
mol_id.py – Molecule visualizer with online InChIKey support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author  : Thomas J. Burton – Savoie Research Group, UND
Updated : 2025-07-22
License : MIT

Description
-----------
Accepts a SMILES string, InChI string, or InChIKey and generates a canonical
SMILES and 2D molecule image. SMILES and InChI are resolved locally using RDKit.
InChIKeys are resolved online via PubChem.

Usage
-----
python mol_id.py "CCO"
python mol_id.py "InChI=1S/CH4/h1H4"
python mol_id.py "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

Options
-------
--no-multicolor   Use black & white atoms
-o                Custom output file name

Dependencies
------------
conda install -c conda-forge rdkit
pip install requests
"""

import argparse
import re
import requests
import textwrap
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DCairo

__version__ = "0.1.0"

# ------------------------------------------------------------------------------
# ACS 1996 Element Color Palette
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
# Format Detection
# ------------------------------------------------------------------------------

def identify_input_type(query: str) -> str:
    """Detect if input is SMILES, InChI, or InChIKey."""
    if query.startswith("InChI="):
        return "inchi"
    elif re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", query):
        return "inchikey"
    else:
        return "smiles"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

def sanitize_filename(text: str) -> str:
    """Make safe filename from SMILES/InChI."""
    return re.sub(r'[^a-zA-Z0-9]', '_', text)

from rdkit.Chem.Draw import rdMolDraw2D

def _apply_palette(opts, mol, multicolor: bool):
    """Try RDKit palette APIs in a version-agnostic way (no halos)."""
    if not multicolor:
        # crisp B/W
        if hasattr(opts, "useBWAtomPalette"):
            opts.useBWAtomPalette()
        return

    # Build an element (Z -> color) palette from your ACS_1996_COLORS
    elem_palette = {}
    for a in mol.GetAtoms():
        col = ACS_1996_COLORS.get(a.GetSymbol())
        if col is not None:
            elem_palette[a.GetAtomicNum()] = col

    # Try the known RDKit attribute names across versions
    if hasattr(opts, "colourPalette"):
        # Newer builds
        for z, col in elem_palette.items():
            opts.colourPalette[z] = col
        return
    if hasattr(opts, "atomColourPalette"):
        # Some builds use this name
        for z, col in elem_palette.items():
            opts.atomColourPalette[z] = col
        return
    if hasattr(opts, "atomColours"):
        # Older/variant builds (per-atom dict)
        # Fallback if neither element-palette exists:
        atom_palette = {a.GetIdx(): elem_palette.get(a.GetAtomicNum(), (0,0,0)) for a in mol.GetAtoms()}
        try:
            opts.atomColours = atom_palette
        except Exception:
            # If even that fails, just leave defaults (still looks fine)
            pass

def draw_molecule(smiles: str, output_file: str, multicolor: bool = True):
    """Draw and save a 2D depiction of a molecule (no highlight halos)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    rdDepictor.Compute2DCoords(mol)

    drawer = MolDraw2DCairo(300, 300)
    opts = drawer.drawOptions()

    # Clean, publication-ish look
    try:
        opts.setBackgroundColour((1, 1, 1, 1))
    except Exception:
        pass
    for attr, val in (("bondLineWidth", 2), ("fixedBondLength", 25.0)):
        if hasattr(opts, attr):
            try:
                setattr(opts, attr, val)
            except Exception:
                pass

    # Apply palette compatibly (prefers per-element, no halos)
    _apply_palette(opts, mol, multicolor)

    # Standard draw path (no highlighting)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    with open(output_file, "wb") as f:
        f.write(drawer.GetDrawingText())
    print(f"[✓] Image saved to: {output_file}")
    
    if args.show:
        try:
            from IPython import get_ipython
            ip = get_ipython()
            if ip is not None:
                from IPython.display import Image as _IPyImage, display as _ipydisplay
                _ipydisplay(_IPyImage(filename=output_path))
        except Exception:
            # Fail silently – keep CLI behavior unchanged
            pass


def resolve_inchikey_to_smiles(inchikey: str) -> tuple[str, int, str] | None:
    """
    Resolve InChIKey to SMILES using CID → PubChem SMILES → RDKit canonical SMILES.

    Returns
    -------
    tuple
        (rdkit_canonical_smiles, cid, pubchem_smiles)
    """
    try:
        # Step 1: resolve CID from InChIKey
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
        cid_resp = requests.get(cid_url, timeout=10)
        cid_resp.raise_for_status()
        cid = cid_resp.json()["IdentifierList"]["CID"][0]

        # Step 2: get PubChem Canonical SMILES
        smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
        smiles_resp = requests.get(smiles_url, timeout=10)
        smiles_resp.raise_for_status()

        props = smiles_resp.json()["PropertyTable"]["Properties"][0]
        pubchem_smiles = props["CanonicalSMILES"]

        # Step 3: Canonicalize via RDKit
        mol = Chem.MolFromSmiles(pubchem_smiles)
        if mol:
            rdkit_smiles = Chem.MolToSmiles(mol, canonical=True)
            return rdkit_smiles, cid, pubchem_smiles
        else:
            print(f"[!] PubChem SMILES could not be parsed by RDKit: {pubchem_smiles}")
            return None

    except Exception as e:
        print(f"[!] PubChem fallback failed: {e}")
        return None



def get_smiles_from_input(input_str: str) -> tuple[str, str, int | None]:
    """
    Get RDKit-canonical SMILES from input (SMILES, InChI, or InChIKey),
    and return original as well.

    Returns
    -------
    tuple
        (rdkit_smiles, original_input_smiles, cid or None)
    """
    input_type = identify_input_type(input_str)

    if input_type == "smiles":
        mol = Chem.MolFromSmiles(input_str)
        if not mol:
            raise ValueError(f"Invalid SMILES: {input_str}")
        rdkit_smiles = Chem.MolToSmiles(mol, canonical=True)
        return rdkit_smiles, input_str, resolve_smiles_or_inchi_to_cid(input_str, "smiles")

    elif input_type == "inchi":
        mol = Chem.MolFromInchi(input_str)
        if not mol:
            raise ValueError(f"Invalid InChI: {input_str}")
        rdkit_smiles = Chem.MolToSmiles(mol, canonical=True)
        return rdkit_smiles, input_str, resolve_smiles_or_inchi_to_cid(input_str, "inchi")

    elif input_type == "inchikey":
        result = resolve_inchikey_to_smiles(input_str)
        if result:
            rdkit_smiles, cid, pubchem_smiles = result
            return rdkit_smiles, pubchem_smiles, cid
        else:
            raise ValueError(f"Could not resolve InChIKey: {input_str}")

    else:
        raise ValueError("Unsupported input format.")


def get_iupac_name_from_cid(cid: int) -> str | None:
    """
    Retrieve IUPAC name from PubChem CID.
    
    Parameters
    ----------
    cid : int
        PubChem Compound ID.
    
    Returns
    -------
    str or None
        IUPAC name if available.
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json()["PropertyTable"]["Properties"][0]["IUPACName"]
    except Exception as e:
        print(f"[!] IUPAC name lookup failed: {e}")
        return None
    
def get_common_names_from_cid(cid: int, top_n: int = 3) -> list[str]:
    """
    Retrieve the top N most human-friendly common names from PubChem.

    Parameters
    ----------
    cid : int
        PubChem Compound ID.
    top_n : int
        How many names to return.

    Returns
    -------
    list of str
        List of common/trivial names (e.g., glucose, EtOH).
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        response = requests.get(url, timeout=10)

        if response.status_code == 404:
            print("[!] No synonyms available for this CID.")
            return []

        response.raise_for_status()
        synonyms = response.json()["InformationList"]["Information"][0]["Synonym"]

        registry_terms = (
            "CID", "CHEBI", "EINECS", "UNII", "InChI", "IUPAC", "PUBCHEM",
            "CHEMBL", "NSC", "ZINC", "EC ", "BRN", "CAS-"
        )

        # Filter out registry/systematic names
        filtered = [
            s for s in synonyms
            if not any(tag in s.upper() for tag in registry_terms)
            and 2 <= len(s) <= 40
        ]

        # Score and sort
        ranked = sorted(filtered, key=lambda x: (-score_common_name(x), len(x)))
        return ranked[:top_n]

    except Exception as e:
        print(f"[!] Common name lookup failed: {e}")
        return []


def resolve_smiles_or_inchi_to_cid(identifier: str, input_type: str) -> int | None:
    """
    Use PubChem to resolve a SMILES or InChI string to a CID.

    Parameters
    ----------
    identifier : str
        The SMILES or InChI string.
    input_type : str
        Either "smiles" or "inchi".

    Returns
    -------
    int or None
        CID if found, else None.
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{input_type}/{requests.utils.quote(identifier)}/cids/JSON"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json()["IdentifierList"]["CID"][0]
    except Exception as e:
        print(f"[!] CID lookup failed for {input_type}: {e}")
        return None

def score_common_name(name: str) -> int:
    """
    Assign a heuristic score to a synonym to prioritize common/trivial names.

    Parameters
    ----------
    name : str
        A synonym from PubChem.

    Returns
    -------
    int
        Higher scores indicate more desirable (common/trivial) names.
    """
    score = 0
    if name.islower():
        score += 3
    elif name[0].isupper():
        score += 2
    if len(name) <= 8:
        score += 3
    elif len(name) <= 15:
        score += 1
    if not any(c in name for c in "[]/\\()"):
        score += 2
    return score

def get_inchi_and_inchikey_from_cid(cid: int) -> tuple[str, str] | None:
    """
    Retrieve InChI and InChIKey from a PubChem CID.

    Parameters
    ----------
    cid : int
        PubChem Compound ID.

    Returns
    -------
    tuple of (InChI, InChIKey), or None
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChI,InChIKey/JSON"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        props = response.json()["PropertyTable"]["Properties"][0]
        return props["InChI"], props["InChIKey"]
    except Exception as e:
        print(f"[!] InChI/InChIKey lookup failed: {e}")
        return None


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mol_id.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Draw a 2D molecule image from a SMILES, InChI, or InChIKey.",
        epilog=textwrap.dedent(
            """\
            Usage:
            python mol_id.py "CCO"
            python mol_id.py "InChI=1S/CH4/h1H4"
            python mol_id.py "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

            Options:
            --no-multicolor        Use black & white atoms
            -o OUTPUT, --output OUTPUT
                                    Custom output image path

            Examples:
            # Ethanol
            python mol_id.py "CCO"

            # Methane by InChI
            python mol_id.py "InChI=1S/CH4/h1H4"

            # Glucose by InChIKey
            python mol_id.py "WQZGKKKJIJFFOK-GASJEMHNSA-N"
            """
        ),
    )
    p.add_argument("input", help="Molecule string (SMILES, InChI, or InChIKey)")
    p.add_argument("-o", "--output", help="Output image filename (default: <smiles>_viz.png)")
    p.add_argument("--no-multicolor", action="store_true", help="Disable ACS coloring (use black and white)")
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    p.add_argument("--show", action="store_true",
               help="If running inside IPython (e.g. via %run), display the image inline")

    return p


def parse_args():
    return build_parser().parse_args()

# ------------------------------------------------------------------------------
# Entrypoint
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    try:
        input_type = identify_input_type(args.input)

        rdkit_smiles = None
        original_smiles = None
        cid = None

        # --- Handle all input types uniformly ---
        rdkit_smiles, original_smiles, cid = get_smiles_from_input(args.input)

        # --- Show PubChem SMILES if input was InChIKey ---
        if input_type == "inchikey":
            print(f"[✓] PubChem SMILES: {original_smiles}")
            print(f"[✓] RDKit Canonical SMILES: {rdkit_smiles}")

        # --- Show both if original SMILES differs ---
        elif input_type == "smiles" and original_smiles != rdkit_smiles:
            print(f"[✓] Input SMILES (original): {original_smiles}")
            print(f"[✓] RDKit Canonical SMILES: {rdkit_smiles}")

        # --- Show canonicalized SMILES directly ---
        else:
            print(f"[✓] SMILES (canonicalized): {rdkit_smiles}")

        # --- If we have CID, show identifiers and names ---
        if cid:
            inchi_info = get_inchi_and_inchikey_from_cid(cid)
            if inchi_info:
                inchi, inchikey = inchi_info
                print(f"[✓] InChI: {inchi}")
                print(f"[✓] InChIKey: {inchikey}")

            iupac = get_iupac_name_from_cid(cid)
            if iupac:
                print(f"[✓] IUPAC name: {iupac}")

            common_names = get_common_names_from_cid(cid)
            if common_names:
                print(f"[✓] Common names: {', '.join(common_names)}")

        # --- Draw image from RDKit canonical SMILES ---
        output_path = args.output or f"{sanitize_filename(rdkit_smiles)}_viz.png"
        draw_molecule(rdkit_smiles, output_path, multicolor=not args.no_multicolor)

    except Exception as e:
        print(f"[✗] Error: {e}")
