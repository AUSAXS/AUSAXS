#!/usr/bin/env python3
"""Build the figures shown on the crystals and powder diffraction application pages.

Calculates the scattering of the silver nanoparticle Ag_crystal.cif with pyAUSAXS, exactly as in the pages' examples,
and draws it in the site's colours. Writes two pairs to docs/src/assets/figures/:

    crystal-debye.svg, crystal-debye-dark.svg    I(q) up to q = 10 Å⁻¹ (crystals.md)
    powder-pattern.svg, powder-pattern-dark.svg  I(2θ) for Cu Kα radiation (powder-diffraction.md)

Use each as a pair:

    ![alt](../assets/figures/<name>.svg#only-light)
    ![alt](../assets/figures/<name>-dark.svg#only-dark)

Needs pyausaxs and matplotlib. Usage: python docs/scripts/build_debye_figures.py
"""
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import numpy as np
import pyausaxs as ausaxs

from build_figures import OUTPUT
from build_fit_figure import THEMES

DOCS = Path(__file__).resolve().parent.parent
DATA = DOCS.parent / "tests" / "files"
WAVELENGTH = 1.5406  # Cu Kα, in Å


def calculate() -> dict:
    ausaxs.settings.histogram(bin_width=0.05) # the default of ausaxs debye
    mol = ausaxs.create_molecule(str(DATA / "Ag_crystal.cif"))
    q = np.linspace(0.01, 10, 2000)
    _, I = mol.debye_raw(q, form_factors=True)

    q_powder = np.linspace(0.01, 4*np.pi/WAVELENGTH, 2000)
    _, I_powder = mol.debye_raw(q_powder, form_factors=True)
    two_theta = np.degrees(2*np.arcsin(q_powder*WAVELENGTH/(4*np.pi)))
    return {"q": q, "I": I, "two_theta": two_theta, "I_powder": I_powder}


def style(colours: dict) -> None:
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


def axes(colours: dict):
    style(colours)
    fig, ax = plt.subplots(figsize=(7, 3.4))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.15)
    fig.patch.set_alpha(0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.patch.set_alpha(0)
    ax.grid(True, color=colours["line"], lw=0.6)
    ax.set_axisbelow(True)
    return fig, ax


def draw_crystal(result: dict, colours: dict):
    fig, ax = axes(colours)
    ax.plot(result["q"], result["I"], color=colours["fit"], lw=1.6)
    ax.set_yscale("log")
    ax.set_xlim(0, 10)
    ax.set_xlabel(r"$q$ (Å$^{-1}$)")
    ax.set_ylabel("$I(q)$")
    return fig


def draw_powder(result: dict, colours: dict):
    fig, ax = axes(colours)
    window = result["two_theta"] >= 25
    ax.plot(result["two_theta"][window], result["I_powder"][window], color=colours["fit"], lw=1.6)
    ax.set_xlim(25, 150)
    ax.set_ylim(0, None)
    ax.set_yticklabels([])
    ax.set_xlabel(r"$2\theta$ (°), Cu K$\alpha$")
    ax.set_ylabel("$I$")
    return fig


def main() -> int:
    result = calculate()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, draw in [("crystal-debye", draw_crystal), ("powder-pattern", draw_powder)]:
        for scheme, colours in THEMES.items():
            fig = draw(result, colours)
            buffer = io.StringIO()
            fig.savefig(buffer, format="svg", metadata={"Date": None})
            plt.close(fig)
            out = OUTPUT / (f"{name}.svg" if scheme == "light" else f"{name}-dark.svg")
            out.write_text(buffer.getvalue(), encoding="utf-8")
            print(f"-> {out.relative_to(DOCS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
