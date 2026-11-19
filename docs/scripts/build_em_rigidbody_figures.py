#!/usr/bin/env python3
"""Build the figures shown on the EM validation and rigid-body refinement application pages.

Runs the same command-line calls as the pages' examples, on the SARS-CoV-2 main protease (SASBDB SASDJG5), and
draws the results in the site's colours. Writes two pairs to docs/src/assets/figures/:

    em-validation.svg, em-validation-dark.svg  the chi2 landscape and the fit at its minimum (em.md)
    rigidbody-fit.svg, rigidbody-fit-dark.svg  the fit before and after refinement (rigid-body-refinement.md)

Use each as a pair:

    ![alt](../assets/figures/<name>.svg#only-light)
    ![alt](../assets/figures/<name>-dark.svg#only-dark)

The data are read from tests/files, which the pages link to, and the configuration script from the code block on
rigid-body-refinement.md. The rigid-body refinement is stochastic, so its figure changes slightly with every build.

Needs the ausaxs command-line tool and matplotlib. Set AUSAXS to the executable if it is not on the path.
Usage: python docs/scripts/build_em_rigidbody_figures.py
"""
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import numpy as np

from build_figures import OUTPUT
from build_fit_figure import THEMES

DOCS = Path(__file__).resolve().parent.parent
DATA = DOCS.parent / "tests" / "files"
MAP = DATA / "emd_24889.map"
RIGIDBODY_PAGE = DOCS / "src" / "rigid-body-refinement.md"
AUSAXS = os.environ.get("AUSAXS", "ausaxs")


def read_blocks(path: Path) -> list[np.ndarray]:
    """The numeric rows of a .fit or .plot file, split into one array per block of consecutive rows."""
    blocks, rows = [], []
    for line in path.read_text().splitlines() + [""]:
        try:
            values = [float(x) for x in line.split()]
        except ValueError:
            values = []
        if values:
            rows.append(values)
        elif rows:
            blocks.append(np.array(rows))
            rows = []
    return blocks


def chi2r(path: Path) -> float:
    return float(next(l for l in path.read_text().splitlines() if l.startswith("chi2=")).split()[0][5:])


def run(args: list[str], cwd: Path) -> None:
    subprocess.run([AUSAXS, *args], cwd=cwd, check=True, stdout=subprocess.DEVNULL)


def run_em(work: Path) -> dict:
    shutil.copy(DATA / "SASDJG5.dat", work)
    shutil.copy(MAP, work)
    run(["em", "emd_24889.map", "SASDJG5.dat"], work)
    out = work / "output" / "em" / "SASDJG5" / "emd_24889"
    smooth, points, minimum = read_blocks(out / "chi2_evaluated_points_limited.png.plot")[:3]
    return {"smooth": smooth, "points": points, "minimum": minimum,
            "fit": read_blocks(out / "ausaxs.fit")[0], "chi2r": chi2r(out / "ausaxs.fit")}


def run_rigidbody(work: Path) -> dict:
    for name in ("SASDJG5.dat", "SASDJG5_single.pdb"):
        shutil.copy(DATA / name, work)
    # the script shown on the page, so the figure always matches it
    script = re.search(r"^```\n(output output/dimer\n.*?)^```", RIGIDBODY_PAGE.read_text(), re.S | re.M).group(1)
    (work / "dimer.conf").write_text(script)
    run(["rigidbody", "dimer.conf"], work)
    out = work / "output" / "dimer"
    return {name: read_blocks(out / f"{name}.fit")[0] for name in ("initial", "final")} | \
           {f"{name}_chi2r": chi2r(out / f"{name}.fit") for name in ("initial", "final")}


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


def finish(ax, colours: dict) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.patch.set_alpha(0)
    ax.grid(True, color=colours["line"], lw=0.6)
    ax.set_axisbelow(True)


def legend(ax, colours: dict, reverse: bool = True) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if reverse:
        handles, labels = handles[::-1], labels[::-1]
    for text in ax.legend(handles, labels, frameon=False, loc="upper right", handlelength=1.6).get_texts():
        text.set_color(colours["ink"])


def data_style(colours: dict) -> dict:
    return dict(fmt="o", ms=2.5, color=colours["data"], ecolor=colours["data"], elinewidth=0.6, alpha=0.6)


def draw_em(result: dict, colours: dict):
    style(colours)
    fig, (left, right) = plt.subplots(1, 2, figsize=(7, 3.2), gridspec_kw={"wspace": 0.3})
    fig.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.16)
    fig.patch.set_alpha(0)

    points, smooth, minimum = result["points"], result["smooth"], result["minimum"]
    left.plot(points[:, 0], points[:, 1], "o", ms=3, color=colours["data"], alpha=0.75, label="evaluated levels")
    left.plot(smooth[:, 0], smooth[:, 1], color=colours["fit"], lw=1.6, label="rolling average")
    left.axvline(minimum[0, 0], color=colours["fit"], lw=1, ls="--")
    left.set_ylim(0, None)
    left.set_xlabel(r"threshold level ($\sigma$)")
    left.set_ylabel(r"$\chi^2_r$")
    legend(left, colours, reverse=False)

    q, I, Ierr, model = (result["fit"][:, i] for i in range(4))
    right.errorbar(q, I, yerr=Ierr, label="SASDJG5.dat", **data_style(colours))
    right.plot(q, model, color=colours["fit"], lw=1.8, zorder=3, label=rf"EMD-24889, $\chi^2_r = {result['chi2r']:.1f}$")
    right.set_yscale("log")
    right.set_ylim(0.3*model.min(), None) # a few noisy points near zero would otherwise stretch the axis
    right.set_xlabel(r"$q$ (Å$^{-1}$)")
    right.set_ylabel("$I(q)$")
    legend(right, colours)

    for ax in (left, right):
        finish(ax, colours)
    return fig


def draw_rigidbody(result: dict, colours: dict):
    style(colours)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7, 4.6), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})
    fig.subplots_adjust(left=0.1, right=0.98, top=0.97, bottom=0.11)
    fig.patch.set_alpha(0)

    initial, final = result["initial"], result["final"]
    q, I, Ierr = final[:, 0], final[:, 1], final[:, 2]
    top.errorbar(q, I, yerr=Ierr, label="SASDJG5.dat", **data_style(colours))
    top.plot(q, initial[:, 3], color=colours["muted"], lw=1.4, ls="--", label=rf"initial, $\chi^2_r = {result['initial_chi2r']:.0f}$")
    top.plot(q, final[:, 3], color=colours["fit"], lw=2, zorder=3, label=rf"refined, $\chi^2_r = {result['final_chi2r']:.1f}$")
    top.set_yscale("log")
    top.set_ylim(0.3*min(initial[:, 3].min(), final[:, 3].min()), None) # a few noisy points near zero would otherwise stretch the axis
    top.set_ylabel("$I(q)$")
    legend(top, colours)

    residuals = (I - final[:, 3])/Ierr
    bottom.axhline(0, color=colours["fit"], lw=1.2)
    bottom.plot(q, residuals, "o", ms=2.5, color=colours["data"], alpha=0.6)
    lim = max(4, 1.15*np.abs(residuals).max())
    bottom.set_ylim(-lim, lim)
    bottom.set_ylabel(r"$\Delta I / \sigma$")
    bottom.set_xlabel(r"$q$ (Å$^{-1}$)")
    bottom.set_xlim(0, q.max()*1.02)

    for ax in (top, bottom):
        finish(ax, colours)
    return fig


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "em").mkdir()
        (Path(tmp) / "rigidbody").mkdir()
        results = {"em-validation": (run_em(Path(tmp) / "em"), draw_em),
                   "rigidbody-fit": (run_rigidbody(Path(tmp) / "rigidbody"), draw_rigidbody)}

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, (result, draw) in results.items():
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
