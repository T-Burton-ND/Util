#!/usr/bin/env python3
"""
color_grad.py

Generate a smooth color gradient between two hex codes and optionally save a PNG
showing swatches with hex labels.

Usage:
  python color_grad.py "#071236" "#5FECFF" 10 --png
  python color_grad.py "#071236" "#5FECFF" 10 --space rgb --png --out my_palette.png
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ---------- sRGB <-> OKLab (Björn Ottosson) ----------
# sRGB in 0..1

def _srgb_to_linear(c):
    c = np.asarray(c)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)

def _linear_to_srgb(c):
    c = np.asarray(c)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * (np.clip(c, 0, None) ** (1/2.4)) - 0.055)

def _rgb_to_oklab(rgb):
    # rgb: (..., 3) in 0..1
    r, g, b = _srgb_to_linear(rgb[..., 0]), _srgb_to_linear(rgb[..., 1]), _srgb_to_linear(rgb[..., 2])

    l_ = 0.4122214708*r + 0.5363325363*g + 0.0514459929*b
    m_ = 0.2119034982*r + 0.6806995451*g + 0.1073969566*b
    s_ = 0.0883024619*r + 0.2817188376*g + 0.6299787005*b

    l_c = np.cbrt(l_)
    m_c = np.cbrt(m_)
    s_c = np.cbrt(s_)

    L = 0.2104542553*l_c + 0.7936177850*m_c - 0.0040720468*s_c
    a = 1.9779984951*l_c - 2.4285922050*m_c + 0.4505937099*s_c
    b = 0.0259040371*l_c + 0.7827717662*m_c - 0.8086757660*s_c

    return np.stack([L, a, b], axis=-1)

def _oklab_to_rgb(lab):
    L, a, b = lab[..., 0], lab[..., 1], lab[..., 2]

    l_c = L + 0.3963377774*a + 0.2158037573*b
    m_c = L - 0.1055613458*a - 0.0638541728*b
    s_c = L - 0.0894841775*a - 1.2914855480*b

    l_ = l_c**3
    m_ = m_c**3
    s_ = s_c**3

    r = +4.0767416621*l_ - 3.3077115913*m_ + 0.2309699292*s_
    g = -1.2684380046*l_ + 2.6097574011*m_ - 0.3413193965*s_
    b = -0.0041960863*l_ - 0.7034186147*m_ + 1.7076147010*s_

    rgb = np.stack([r, g, b], axis=-1)
    return np.clip(_linear_to_srgb(rgb), 0.0, 1.0)

# ---------- Utilities ----------

def hex_to_rgb(hex_color):
    return np.array(mcolors.to_rgb(hex_color))

def rgb_to_hex(rgb):
    return mcolors.to_hex(np.clip(rgb, 0, 1), keep_alpha=False)

def gradient_colors(start_hex, end_hex, n=10, space="oklab"):
    start_rgb = hex_to_rgb(start_hex)
    end_rgb = hex_to_rgb(end_hex)

    if space.lower() == "oklab":
        a = _rgb_to_oklab(start_rgb[None, :])[0]
        b = _rgb_to_oklab(end_rgb[None, :])[0]
        labs = np.linspace(a, b, n)
        rgbs = _oklab_to_rgb(labs)
    else:
        rgbs = np.linspace(start_rgb, end_rgb, n)

    return [rgb_to_hex(c) for c in rgbs]

def relative_luminance(rgb):
    # WCAG-ish luminance for sRGB 0..1
    def lin(c):
        return np.where(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    r, g, b = lin(np.array(rgb))
    return 0.2126*r + 0.7152*g + 0.0722*b

def save_palette_png(colors, path, cell_w=240, cell_h=120, dpi=200):
    """
    Draw a vertical strip of color cells with hex labels and save as PNG.
    """
    n = len(colors)
    fig_w, fig_h = cell_w / dpi, (cell_h * n) / dpi
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n)
    ax.axis("off")

    for i, hexc in enumerate(colors):
        ax.add_patch(plt.Rectangle((0, i), 1, 1, color=hexc, ec="black", lw=0.5))
        rgb = np.array(mcolors.to_rgb(hexc))
        text_color = "white" if relative_luminance(rgb) < 0.5 else "black"
        ax.text(
            0.5, i + 0.5, hexc.upper(),
            ha="center", va="center",
            fontsize=cell_h * 0.25 / (dpi / 72.0),
            color=text_color, fontweight="bold"
        )

    plt.tight_layout(pad=0)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def auto_name_png(start_hex, end_hex, n):
    return f"palette_{start_hex.lstrip('#')}_to_{end_hex.lstrip('#')}_n{n}.png"

# ---------- CLI ----------

def parse_args():
    p = argparse.ArgumentParser(description="Generate a color gradient and optional PNG swatch bar.")
    p.add_argument("start_hex", help='Start color, e.g. "#071236"')
    p.add_argument("end_hex", help='End color, e.g. "#5FECFF"')
    p.add_argument("n", type=int, help="Number of colors (including endpoints)")
    p.add_argument("--space", choices=["oklab", "rgb"], default="oklab",
                   help="Interpolation space (default: oklab)")
    p.add_argument("--png", action="store_true", help="Save a PNG swatch bar", default=True)
    p.add_argument("--out", default=None, help="Output PNG filename (used if --png)")
    p.add_argument("--txt", default=None, help="Also write hex list to a text file")
    p.add_argument("--cell-w", type=int, default=240, help="Cell width in px")
    p.add_argument("--cell-h", type=int, default=120, help="Cell height in px")
    p.add_argument("--dpi", type=int, default=200, help="PNG DPI")
    return p.parse_args()

def main():
    args = parse_args()

    colors = gradient_colors(args.start_hex, args.end_hex, n=args.n, space=args.space)

    print("\nGenerated colors:")
    for c in colors:
        print(c)

    if args.txt:
        with open(args.txt, "w") as f:
            for c in colors:
                f.write(c + "\n")
        print(f"\nWrote hex list to: {args.txt}")

    if args.png:
        out = args.out or auto_name_png(args.start_hex, args.end_hex, args.n)
        save_palette_png(colors, out, cell_w=args.cell_w, cell_h=args.cell_h, dpi=args.dpi)
        print(f"Wrote swatch PNG to: {out}")

if __name__ == "__main__":
    main()
