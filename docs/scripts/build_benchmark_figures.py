#!/usr/bin/env python3
"""Build the benchmark figures from the results of tests/benchmark/benchmark_paper.

Reads results.jsonl and iq/ from the benchmark's output folder (default output/paper_benchmark) and draws three
figures in the site's colours. Writes three pairs to docs/src/assets/figures/:

    benchmark-scaling.svg, benchmark-scaling-dark.svg      time against size: SAXS fits and the pure Debye sum
    benchmark-scaling-sas.svg, ...-dark.svg                 the SAXS fits only
    benchmark-scaling-debye.svg, ...-dark.svg               the pure Debye sum only
    benchmark-accuracy.svg, benchmark-accuracy-dark.svg    deviation of the Debye sum from the exact sum
    benchmark-gradient.svg, benchmark-gradient-dark.svg    time of the gradient of the Debye sum against size
    benchmark-scaling-modern.svg, ...-dark.svg              AUSAXS alone on a modern desktop, from LOCAL_RESULTS

Use each as a pair:

    ![alt](assets/figures/<name>.svg#only-light)
    ![alt](assets/figures/<name>-dark.svg#only-dark)

Points are the mean of the timed runs. A single run over the benchmark's time budget ends its series; it is drawn as
an open marker without a line. Deviations are |I/I_exact - 1| after removing each program's constant normalisation
(AUSAXS sums over ordered pairs, twice DebyeCalculator's sum); the reference is DebyeCalculator's direct sum in float64.

Needs numpy and matplotlib. Usage: python docs/scripts/build_benchmark_figures.py [output/paper_benchmark [local results]]
"""
import io
import json
import statistics
import sys
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.lines import Line2D

from build_figures import OUTPUT

DOCS = Path(__file__).resolve().parent.parent
RESULTS = DOCS.parent / "output" / "paper_benchmark"
# The AUSAXS-only [debye] run on a modern desktop (AMD Ryzen 5 9600X, AMD Radeon RX 9070 XT), with competitors off and
# AUSAXS_BENCH_CPU_MAX_ATOMS=300000. Its last GPU sizes are missing: the run took down the desktop the GPU also drives.
LOCAL_RESULTS = DOCS.parent / "output" / "paper_benchmark_local"

# Site colours from docs/src/stylesheets/extra.css; AUSAXS wears the site's gold. Programs that are compared directly
# (the three SAS competitors; the two Debye competitors) are separated under protan, deutan and tritan simulation (OKLab
# ΔE ≥ 9.4 against each other and AUSAXS). Every program also has its own marker, so identity never rests on colour.
THEMES = {
    "light": {"ink": "#1f2733", "muted": "#57606a", "line": "#d3dae3",
              "ausaxs": "#b07c12", "pepsi-saxs": "#2a78d6", "foxs": "#e87ba4", "crysol": "#3a4556",
              "debyecalculator": "#6250d6", "aes-debye": "#1baf7a"},
    "dark":  {"ink": "#e6edf3", "muted": "#8b949e", "line": "#3a4556",
              "ausaxs": "#e3b341", "pepsi-saxs": "#3987e5", "foxs": "#d55181", "crysol": "#c9d1d9",
              "debyecalculator": "#9085e9", "aes-debye": "#199e70"},
}
MARKER = {"ausaxs": "o", "pepsi-saxs": "s", "foxs": "^", "crysol": "D", "debyecalculator": "v", "aes-debye": "P"}

# Line style by kind of calculation: (line style, filled marker)
FIT, CPU, GPU = ("-", True), ("--", True), (":", False)

# The AUSAXS Debye sum is timed at 0.25 and 0.05 Å. 0.05 Å is shown, the default for crystals and powder diffraction.
DEBYE_TIMING = "forward-cpu-bw0.05", "forward-gpu-bw0.05"
# The AUSAXS accuracy at the default bin width for crystals and powder diffraction.
ACCURACY_AUSAXS = "bw0.05"
ACCURACY_SIZE = "Ag_R30"


def load(folder: Path) -> list[dict]:
    records = {}
    with open(folder / "results.jsonl") as f:
        for line in f:
            r = json.loads(line)
            records[r["key"]] = r  # the last record of a key wins, so a re-measurement replaces the old one
    return list(records.values())


def series(records, suite, program, variant):
    """(atoms, mean time) of the completed points, and the point over budget that ended the series, if any."""
    match = [r for r in records if (r["suite"], r["program"], r["variant"]) == (suite, program, variant)]
    pts = sorted((r["atoms"], statistics.fmean(r["samples"])) for r in match if r["status"] == "ok")
    over = [(r["atoms"], r["samples"][0]) for r in match if r["status"] == "over budget"]
    return np.array(pts).T, (over[0] if over else None)


def deviation(folder: Path, program: str, variant: str, structure: str):
    """q and |I/I_exact - 1|, with the program's normalisation convention (1/2, 1 or 2) removed."""
    iq = np.loadtxt(folder / "iq" / f"{program}_{variant}_{structure}.dat")
    ref = np.loadtxt(folder / "iq" / f"debyecalculator_exact-float64_{structure}.dat")
    if iq.shape != ref.shape or np.max(np.abs(iq[:, 0] - ref[:, 0])) > 1e-5:
        raise ValueError(f"{program} {variant} {structure}: q grids differ")
    ratio = iq[:, 1]/ref[:, 1]
    conv = min((0.5, 1.0, 2.0), key=lambda c: abs(ratio[0]/c - 1))
    if abs(ratio[0]/conv - 1) > 0.01:
        raise ValueError(f"{program} {variant} {structure}: no known normalisation, I/I_exact = {ratio[0]:.4f}")
    return iq[:, 0], np.abs(ratio/conv - 1)


def rolling_median(y: np.ndarray, window: int) -> np.ndarray:
    half = window//2
    padded = np.pad(y, half, mode="edge")
    return np.median(np.lib.stride_tricks.sliding_window_view(padded, window), axis=1)


def figure(c: dict, size=(8, 3.9)):
    plt.rcParams.update({
        "font.family": ["Inter", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9.5,
        "svg.fonttype": "none",
        "svg.hashsalt": "ausaxs-docs",
        "text.color": c["ink"],
        "axes.labelcolor": c["ink"],
        "axes.edgecolor": c["muted"],
        "xtick.color": c["muted"],
        "ytick.color": c["muted"],
        "xtick.labelcolor": c["ink"],
        "ytick.labelcolor": c["ink"],
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "legend.handlelength": 2.8,
        "legend.title_fontsize": 8.5,
    })
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_alpha(0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.patch.set_alpha(0)
    ax.grid(True, which="major", color=c["line"], lw=0.6)
    ax.set_axisbelow(True)
    return fig, ax


def plot(ax, c, data, program, label, kind=FIT, lw=1.5, marker=True):
    (x, y), over = data
    col = c[program]
    ls, filled = kind
    m = MARKER[program] if marker else "none"
    ax.plot(x, y, color=col, ls=ls, lw=lw, marker=m, markersize=5.5, markeredgewidth=1.1,
            markerfacecolor=col if filled else "none", label=label)
    if over is not None and marker:
        ax.plot(*over, color=col, marker=m, markersize=5.5, markeredgewidth=1.1, markerfacecolor="none", ls="none")


def log_decades(axis):
    axis.set_major_locator(ticker.LogLocator(base=10, numticks=12))
    axis.set_major_formatter(ticker.LogFormatterMathtext())
    axis.set_minor_locator(ticker.NullLocator())


def title_handle():
    return Line2D([], [], ls="none", marker="none")


def fit_series(ax, c, records):
    """Every program runs as a process from the command line, AUSAXS included."""
    plot(ax, c, series(records, "sas", "ausaxs-cpu", "cli"), "ausaxs", "AUSAXS")
    plot(ax, c, series(records, "sas", "pepsi-saxs", "cli"), "pepsi-saxs", "Pepsi-SAXS")
    plot(ax, c, series(records, "sas", "foxs", "cli"), "foxs", "FoXS")
    plot(ax, c, series(records, "sas", "crysol", "cli"), "crysol", "CRYSOL")


def debye_series(ax, c, records):
    plot(ax, c, series(records, "debye", "ausaxs", DEBYE_TIMING[0]), "ausaxs", "AUSAXS, CPU", CPU)
    plot(ax, c, series(records, "debye", "ausaxs", DEBYE_TIMING[1]), "ausaxs", "AUSAXS, GPU", GPU)
    plot(ax, c, series(records, "debye", "aes-debye", "forward-cpu"), "aes-debye", "AES-Debye, CPU", CPU)
    plot(ax, c, series(records, "debye", "aes-debye", "forward-cuda"), "aes-debye", "AES-Debye, GPU", GPU)
    plot(ax, c, series(records, "debye", "debyecalculator", "forward-cpu"), "debyecalculator", "DebyeCalculator, CPU", CPU)
    plot(ax, c, series(records, "debye", "debyecalculator", "forward-cuda"), "debyecalculator", "DebyeCalculator, GPU", GPU)


def time_axes(ax, xlim, ylim):
    ax.set_xscale("log")
    ax.set_yscale("log")
    log_decades(ax.xaxis)
    log_decades(ax.yaxis)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Number of atoms")
    ax.set_ylabel("Time (s)")


def draw_scaling(records, folder, c):
    """Everything on one axis: SAXS fits and the Debye sum."""
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    fit_series(ax, c, records)
    n = len(ax.get_legend_handles_labels()[0])
    debye_series(ax, c, records)
    handles, labels = ax.get_legend_handles_labels()
    time_axes(ax, (1.2e2, 1.5e6), (1e-4, 3e3))
    first = ax.legend(handles[:n], labels[:n], title="Fitting SAXS data, CPU", alignment="left", loc="upper left")
    ax.add_artist(first)
    ax.legend(handles[n:], labels[n:], title="Debye sum, nanoparticles", alignment="left", loc="upper left",
              bbox_to_anchor=(0.2, 1.0))
    return fig


def draw_scaling_sas(records, folder, c):
    """SAXS fits only."""
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    fit_series(ax, c, records)
    time_axes(ax, (4.5e2, 6e4), (3e-3, 1e2))
    ax.legend(title="Fitting SAXS data, CPU", alignment="left", loc="upper left")
    return fig


def draw_scaling_debye(records, folder, c):
    """The Debye sum only."""
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    debye_series(ax, c, records)
    time_axes(ax, (1.2e2, 1.5e6), (2e-4, 3e3))
    ax.legend(title="Debye sum, silver nanoparticles", alignment="left", loc="upper left")
    return fig


def draw_accuracy(records, folder, c):
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    atoms = {r["structure"]: r["atoms"] for r in records if r["suite"] == "accuracy"}

    # a running median, so the spikes at the minima of I(q) don't hide the level
    window = 25
    for program, variant, label, ls in (
        ("aes-debye", "res1.0", "AES-Debye", "-"),
        ("debyecalculator", "float32-f0.12", "DebyeCalculator", "--"),
        ("ausaxs", ACCURACY_AUSAXS, "AUSAXS", "-"),
    ):
        q, d = deviation(folder, program, variant, ACCURACY_SIZE)
        ax.plot(q, rolling_median(d, window), color=c[program], lw=1.6, ls=ls, label=label)
    handles, labels = ax.get_legend_handles_labels()

    ax.set_yscale("log")
    log_decades(ax.yaxis)
    ax.set_xlim(0, 10)
    ax.set_ylim(1e-8, 1e-3)
    ax.set_xlabel("q (Å⁻¹)")
    ax.set_ylabel("|I / I$_\\mathrm{exact}$ − 1|")
    ax.legend(handles[::-1], labels[::-1], loc="upper right", ncol=3,
              title=f"Silver nanoparticle, {atoms[ACCURACY_SIZE]:,} atoms, default settings", alignment="right")
    return fig


def draw_gradient(records, folder, c):
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)

    ours = series(records, "debye", "ausaxs", "gradient-cpu")
    dc_cpu = series(records, "debye", "debyecalculator", "gradient-cpu")
    dc_gpu = series(records, "debye", "debyecalculator", "gradient-cuda")
    plot(ax, c, ours, "ausaxs", "AUSAXS, CPU")
    plot(ax, c, dc_cpu, "debyecalculator", "DebyeCalculator, CPU")
    plot(ax, c, dc_gpu, "debyecalculator", "DebyeCalculator, GPU", GPU)
    plot(ax, c, series(records, "debye", "ausaxs", DEBYE_TIMING[0]), "ausaxs", "AUSAXS, I(q) only", ("-", True),
         lw=0.9, marker=False)

    # the speed-up against the faster of DebyeCalculator's two devices, at the largest size all three completed
    (xa, ya), _ = ours
    best = {}
    for (x, y), _ in (dc_cpu, dc_gpu):
        for n, t in zip(x, y):
            best[n] = min(best.get(n, np.inf), t)
    n = max(v for v in xa if v in best)
    t_ours, t_dc = ya[list(xa).index(n)], best[n]
    ax.annotate("", xy=(n, t_ours*1.25), xytext=(n, t_dc/1.25),
                arrowprops={"arrowstyle": "<->", "color": c["muted"], "lw": 0.9, "shrinkA": 0, "shrinkB": 0})
    ax.annotate(f"{t_dc/t_ours:.0f}× faster", (n, np.sqrt(t_ours*t_dc)), xytext=(7, 0), textcoords="offset points",
                va="center", fontsize=9, color=c["ink"])

    time_axes(ax, (1.2e2, 1.5e6), (1e-4, 3e3))
    ax.legend(loc="upper left", title="Gradient with respect to the atom positions", alignment="left")
    return fig


def draw_scaling_modern(records, local, c):
    """AUSAXS alone on the modern desktop, with the reference machine's AUSAXS in grey."""
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    plot(ax, c, series(local, "debye", "ausaxs", DEBYE_TIMING[0]), "ausaxs", "CPU, bin width 0.05 Å", CPU)
    plot(ax, c, series(local, "debye", "ausaxs", DEBYE_TIMING[1]), "ausaxs", "GPU, bin width 0.05 Å", GPU)
    (x, y), _ = series(local, "debye", "ausaxs", "forward-gpu-bw0.25")
    ax.plot(x, y, color=c["ausaxs"], ls=":", lw=1.5, marker="s", markersize=5, markeredgewidth=1.1,
            markerfacecolor="none", label="GPU, bin width 0.25 Å")
    modern = ax.get_legend_handles_labels()
    for variant, ls, label in ((DEBYE_TIMING[0], "--", "CPU, bin width 0.05 Å"), (DEBYE_TIMING[1], ":", "GPU, bin width 0.05 Å")):
        (x, y), _ = series(records, "debye", "ausaxs", variant)
        ax.plot(x, y, color=c["muted"], ls=ls, lw=1.0, label=label)
    handles, labels = ax.get_legend_handles_labels()
    n = len(modern[0])
    time_axes(ax, (1.2e2, 2.5e6), (1e-4, 3e2))
    first = ax.legend(handles[:n], labels[:n], title="AUSAXS, Ryzen 5 9600X and Radeon RX 9070 XT", alignment="left",
                      loc="upper left")
    ax.add_artist(first)
    ax.legend(handles[n:], labels[n:], title="AUSAXS, i7-8700 and GTX 1070 Ti", alignment="left", loc="upper left",
              bbox_to_anchor=(0.0, 0.7))
    return fig


def main() -> int:
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS
    local_folder = Path(sys.argv[2]) if len(sys.argv) > 2 else LOCAL_RESULTS
    records = load(folder)
    local = load(local_folder)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, draw, arg in [("benchmark-scaling", draw_scaling, folder), ("benchmark-scaling-sas", draw_scaling_sas, folder),
                            ("benchmark-scaling-debye", draw_scaling_debye, folder),
                            ("benchmark-scaling-modern", draw_scaling_modern, local),
                            ("benchmark-accuracy", draw_accuracy, folder), ("benchmark-gradient", draw_gradient, folder)]:
        for scheme, colours in THEMES.items():
            fig = draw(records, arg, colours)
            buffer = io.StringIO()
            fig.savefig(buffer, format="svg", metadata={"Date": None})
            plt.close(fig)
            out = OUTPUT / (f"{name}.svg" if scheme == "light" else f"{name}-dark.svg")
            out.write_text(buffer.getvalue(), encoding="utf-8")
            print(f"-> {out.relative_to(DOCS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
