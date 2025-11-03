#!/usr/bin/env python3
# pathway_pack.py  (arrow-image version)

import argparse, os, io, csv, json
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageColor, ImageChops

# ---------------------------- RDKit drawing helpers ----------------------------
def drop_atom_maps(m):
    m = Chem.Mol(m)
    for a in m.GetAtoms():
        if a.GetAtomMapNum():
            a.SetAtomMapNum(0)
    return m

def make_mol(smi, keep_maps=False):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    if not keep_maps:
        m = drop_atom_maps(m)
    rdDepictor.Compute2DCoords(m)
    return m

def render_png_raw(mol, size=(900,700), bw=False):
    """Render a single molecule to RGBA image."""
    w,h=size
    drawer=rdMolDraw2D.MolDraw2DCairo(w,h)
    opts=drawer.drawOptions()
    try: opts.setBackgroundColour((1,1,1,0))
    except: pass
    if hasattr(opts,"bondLineWidth"): opts.bondLineWidth=2
    if hasattr(opts,"padding"): opts.padding=0.02
    if bw and hasattr(opts,"useBWAtomPalette"): opts.useBWAtomPalette()
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol)
    drawer.FinishDrawing()
    return Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGBA")

def autocrop_rgba(im, margin=6):
    """Trim transparent border and re-add small margin."""
    alpha=im.split()[-1]
    bbox=alpha.getbbox()
    if not bbox: return im
    l,t,r,b=bbox
    l=max(0,l-margin); t=max(0,t-margin)
    r=min(im.width,r+margin); b=min(im.height,b+margin)
    return im.crop((l,t,r,b))

def render_cropped_tile(mol,bw=False):
    base=render_png_raw(mol,size=(900,700),bw=bw)
    return autocrop_rgba(base,margin=8)

# ---------------------------- Composition helpers ------------------------------
def _load_rgba(p):
    im=Image.open(p)
    return im.convert("RGBA") if im.mode!="RGBA" else im

def _scale_to_height(im,H):
    if im.height==H: return im
    new_w=int(im.width*(H/im.height))
    return im.resize((new_w,H),Image.LANCZOS)

def compose_row_dynamic(images,out_path,
                        arrow_img=None,arrow_scale=1.0,
                        gutter="auto",gutter_frac=0.12,min_gutter=40,pad=40,
                        bg="transparent",arrows=True,arrow_width=6,
                        normalize_height=True):
    ims=[_load_rgba(p) for p in images]
    if normalize_height:
        H=max(im.height for im in ims)
        ims=[_scale_to_height(im,H) for im in ims]
    widths=[im.width for im in ims]; H=max(im.height for im in ims)

    # --- Gutters ---
    if gutter=="auto":
        gaps=[max(min_gutter,int(gutter_frac*min(widths[i],widths[i+1])))
              for i in range(len(ims)-1)]
    else:
        g=int(gutter); gaps=[g]*(len(ims)-1)
    total_w=pad*2+sum(widths)+sum(gaps); total_h=pad*2+H

    bg_color=(255,255,255,0) if bg=="transparent" else ImageColor.getrgb(bg)
    if len(bg_color)==3: bg_color=(*bg_color,255)
    canvas=Image.new("RGBA",(total_w,total_h),bg_color)

    # --- Paste molecules ---
    x=pad; centers=[]
    for i,im in enumerate(ims):
        y=pad+(H-im.height)//2
        canvas.paste(im,(x,y),im)
        centers.append((x+im.width//2,pad+H//2))
        x+=im.width+(gaps[i] if i<len(gaps) else 0)

    # --- Arrows ---
    if arrows and len(ims)>1:
        if arrow_img and os.path.exists(arrow_img):
            arrow=Image.open(arrow_img).convert("RGBA")
            w0,h0=arrow.size
            arrow=arrow.resize((int(w0*arrow_scale),int(h0*arrow_scale)),Image.LANCZOS)
            for i in range(len(ims)-1):
                gap=gaps[i]
                gap_w=int(gap*0.8)
                arrow_scaled=arrow.resize(
                    (min(gap_w,arrow.width),arrow.height),Image.LANCZOS)
                left_right=centers[i][0]+widths[i]//2
                right_left=centers[i+1][0]-widths[i+1]//2
                cx=(left_right+right_left)//2
                cy=centers[i][1]
                ax=cx-arrow_scaled.width//2
                ay=cy-arrow_scaled.height//2
                canvas.paste(arrow_scaled,(ax,ay),arrow_scaled)
        else:
            draw=ImageDraw.Draw(canvas)
            for i in range(len(ims)-1):
                left_right=centers[i][0]+widths[i]//2
                right_left=centers[i+1][0]-widths[i+1]//2
                y=centers[i][1]; gap=right_left-left_right
                m=max(10,gap//6)
                start=left_right+m; end=right_left-m
                draw.line([(start,y),(end,y)],fill=(0,0,0,255),width=arrow_width)
                ah=max(10,gap//5); aw=max(6,gap//7)
                head=[(end,y),(end-ah,y-aw),(end-ah,y+aw)]
                draw.polygon(head,fill=(0,0,0,255))

    canvas.save(out_path)
    return out_path

# ----------------------------- Manifest writers -------------------------------
def write_manifest_csv(path,items):
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["index","filename","smiles","width","height"])
        w.writeheader(); [w.writerow(it) for it in items]
def write_manifest_json(path,items):
    with open(path,"w") as f: json.dump(items,f,indent=2)

# ----------------------------------- CLI --------------------------------------
def main():
    ap=argparse.ArgumentParser(description="Render SMILES → cropped PNGs → row image (arrow image optional).")
    ap.add_argument("smiles",nargs="+",help="2–4 SMILES")
    ap.add_argument("-o","--out",default="pathway.png",help="Final PNG")
    ap.add_argument("--outdir",default="mol_pngs",help="Tile directory")
    ap.add_argument("--prefix",default="mol",help="Tile prefix (mol_01.png...)")
    ap.add_argument("--bw",action="store_true",help="Black & white atoms")
    ap.add_argument("--keep-maps",action="store_true",help="Keep atom-map numbers")
    ap.add_argument("--no-arrows",action="store_true",help="Skip arrows")
    ap.add_argument("--arrow-img",help="PNG arrow file to place between molecules", default="arrow.png")
    ap.add_argument("--arrow-scale",type=float,default=1.0,help="Scale factor for arrow image")
    ap.add_argument("--bg",default="transparent",help="Background color")
    ap.add_argument("--no-normalize",action="store_true",help="Keep RDKit-native heights")
    ap.add_argument("--gutter",default="auto",help='"auto" or pixel gap')
    ap.add_argument("--gutter-frac",type=float,default=0.12,help="Gap fraction for auto mode")
    ap.add_argument("--min-gutter",type=int,default=40,help="Minimum gap (px)")
    args=ap.parse_args()

    if not (2<=len(args.smiles)<=4):
        raise SystemExit("Provide between 2 and 4 SMILES.")
    os.makedirs(args.outdir,exist_ok=True)

    manifest=[]; tiles=[]
    for i,smi in enumerate(args.smiles,1):
        mol=make_mol(smi,keep_maps=args.keep_maps)
        tile=render_cropped_tile(mol,bw=args.bw)
        path=os.path.join(args.outdir,f"{args.prefix}_{i:02d}.png")
        tile.save(path)
        manifest.append({"index":i,"filename":os.path.abspath(path),
                         "smiles":smi,"width":tile.width,"height":tile.height})
        tiles.append(path)
    write_manifest_csv(os.path.join(args.outdir,f"{args.prefix}_manifest.csv"),manifest)
    write_manifest_json(os.path.join(args.outdir,f"{args.prefix}_manifest.json"),manifest)

    final=compose_row_dynamic(tiles,args.out,
                              arrow_img=args.arrow_img,
                              arrow_scale=args.arrow_scale,
                              gutter=args.gutter,gutter_frac=args.gutter_frac,
                              min_gutter=args.min_gutter,pad=40,bg=args.bg,
                              arrows=not args.no_arrows,normalize_height=not args.no_normalize)
    print(f"[✓] Final image: {os.path.abspath(final)}")

if __name__=="__main__":
    main()
