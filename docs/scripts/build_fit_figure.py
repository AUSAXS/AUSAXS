#!/usr/bin/env python3
"""Build the lysozyme fit figure shown on the SAXS application page.

Fits 2epe.pdb to 2epe.dat with pyAUSAXS, exactly as in the page's example, and draws the fit in the site's
colours. The stylized molecule from figures/saxs-decomposition.svg is embedded, as vector graphics, in the
bottom-left corner. Writes saxs-fit.svg and saxs-fit-dark.svg to docs/src/assets/figures/; use them as a pair:

    ![alt](assets/figures/saxs-fit.svg#only-light)
    ![alt](assets/figures/saxs-fit-dark.svg#only-dark)

Needs pyausaxs and matplotlib. Usage: python docs/scripts/build_fit_figure.py
"""
import io
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import pyausaxs as ausaxs

from build_figures import OUTPUT, SOURCES, SVG, remove_hidden, remove_unused_defs, strip_editor_data

DOCS = Path(__file__).resolve().parent.parent
DATA = DOCS.parent / "tests" / "files"

# The "Molecule in vacuum" drawing in the decomposition figure: its group id, and the region around it in the
# source's user units (mm). Update the region if the drawing is moved in Inkscape.
MOLECULE_ID = "g170"
MOLECULE_VIEWBOX = (77.0, 16.5, 27.5, 23.5)

# Colours from docs/src/stylesheets/extra.css. Both schemes draw on a transparent background.
THEMES = {
    "light": {"ink": "#1f2733", "muted": "#57606a", "line": "#d3dae3", "data": "#3a4556", "fit": "#b07c12"},
    "dark":  {"ink": "#e6edf3", "muted": "#8b949e", "line": "#3a4556", "data": "#c9d1d9", "fit": "#e3b341"},
}


def run_fit():
    mol = ausaxs.create_molecule(str(DATA / "2epe.pdb"))
    mol.hydrate()
    return mol.fit(str(DATA / "2epe.dat"))


def molecule_svg() -> str:
    """The molecule drawing as a standalone <svg> fragment, cropped to MOLECULE_VIEWBOX, with prefixed ids."""
    root = ET.parse(SOURCES / "saxs-decomposition.svg").getroot()
    target = next(e for e in root.iter() if e.get("id") == MOLECULE_ID)

    # Keep the definitions, the target, and the groups leading to it (their transforms place it).
    parents = {c: p for p in root.iter() for c in p}
    keep = {target}
    e = target
    while e in parents:
        e = parents[e]
        keep.add(e)
    for parent in [root, *(p for p in keep if p is not target)]:
        for child in list(parent):
            if child not in keep and child.tag != f"{{{SVG}}}defs":
                parent.remove(child)

    strip_editor_data(root)
    remove_hidden(root)
    remove_unused_defs(root)
    for attr in ("width", "height", "role", "aria-label"):
        root.attrib.pop(attr, None)
    root.set("viewBox", " ".join(map(str, MOLECULE_VIEWBOX)))
    root.set("preserveAspectRatio", "xMinYMax meet")

    text = ET.tostring(root, encoding="unicode")
    text = re.sub(r'\bid="([^"]+)"', r'id="mol-\1"', text)
    text = re.sub(r"url\(#([^)]+)\)", r"url(#mol-\1)", text)
    return re.sub(r'href="#([^"]+)"', r'href="#mol-\1"', text)


def draw(fit, colours: dict):
    q, I, Ierr = fit.data_curve()
    model = fit.fitted_curve()
    residuals = fit.residuals()

    plt.rcParams.update({
        "font.family": ["Inter", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 10,
        "svg.fonttype": "none",
        "svg.hashsalt": "ausaxs-docs",
        "text.color": colours["ink"],
        "axes.labelcolor": colours["ink"],
        "axes.edgecolor": colours["muted"],
        "xtick.color": colours["muted"],
        "ytick.color": colours["muted"],
        "xtick.labelcolor": colours["ink"],
        "ytick.labelcolor": colours["ink"],
        "axes.linewidth": 0.8,
    })
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7, 4.6), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})
    fig.subplots_adjust(left=0.1, right=0.98, top=0.97, bottom=0.11)
    fig.patch.set_alpha(0)

    data_style = dict(fmt="o", ms=3, color=colours["data"], ecolor=colours["data"], elinewidth=0.8, alpha=0.75)
    top.errorbar(q, I, yerr=Ierr, label="2epe.dat", **data_style)
    top.plot(q, model, color=colours["fit"], lw=2, zorder=3, label=rf"AUSAXS fit, $\chi^2_r = {fit.chi2r():.2f}$")
    # top.set_xscale("log")
    top.set_yscale("log")
    top.set_ylabel("$I(q)$")
    handles, labels = top.get_legend_handles_labels()
    legend = top.legend(handles[::-1], labels[::-1], frameon=False, loc="upper right", handlelength=1.6)
    for text in legend.get_texts():
        text.set_color(colours["ink"])

    bottom.axhline(0, color=colours["fit"], lw=1.2)
    bottom.plot(q, residuals, "o", ms=3, color=colours["data"], alpha=0.75)
    lim = max(4, 1.15 * abs(residuals).max())
    # bottom.set_xscale("log")
    bottom.set_ylim(-lim, lim)
    bottom.set_yticks([-3, 0, 3])
    bottom.set_ylabel(r"$\Delta I / \sigma$")
    bottom.set_xlabel(r"$q$ (Å$^{-1}$)")
    bottom.set_xlim(0, q.max() * 1.02)

    for ax in (top, bottom):
        ax.spines[["top", "right"]].set_visible(False)
        ax.patch.set_alpha(0)
        ax.grid(True, color=colours["line"], lw=0.6)
        ax.set_axisbelow(True)
    return fig, top


def embed_molecule(svg: str, fig, ax, molecule: str) -> str:
    """Place the molecule in the bottom-left corner of ax, in the figure's own SVG coordinates (points)."""
    width_pt, height_pt = fig.get_size_inches() * 72
    box = ax.get_position()
    w = 0.30 * box.width * width_pt
    h = w * MOLECULE_VIEWBOX[3] / MOLECULE_VIEWBOX[2]
    x = box.x0 * width_pt + 0.03 * box.width * width_pt
    y = (1 - box.y0) * height_pt - 0.05 * box.height * height_pt - h
    placed = molecule.replace("<svg ", f'<svg x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" ', 1)
    return svg.replace("</svg>", placed + "\n</svg>")


def main() -> int:
    fit = run_fit()
    print(f"chi2r = {fit.chi2r():.3f}, parameters: {fit.fit_parameters()}")
    molecule = molecule_svg()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for scheme, colours in THEMES.items():
        fig, top = draw(fit, colours)
        buffer = io.StringIO()
        fig.savefig(buffer, format="svg", metadata={"Date": None})
        plt.close(fig)
        out = OUTPUT / ("saxs-fit.svg" if scheme == "light" else "saxs-fit-dark.svg")
        out.write_text(embed_molecule(buffer.getvalue(), fig, top, molecule), encoding="utf-8")
        print(f"-> {out.relative_to(DOCS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
