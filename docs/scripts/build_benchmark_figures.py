#!/usr/bin/env python3
"""Build the benchmark figures from the results of tests/benchmark/benchmark_paper.

Reads results.jsonl and iq/ from the benchmark's output folder (default output/paper_benchmark) and draws the
figures in the site's colours. Writes these pairs to docs/src/assets/figures/:

    benchmark-scaling.svg, benchmark-scaling-dark.svg      time against size: SAXS fits and the pure Debye sum
    benchmark-scaling-sas.svg, ...-dark.svg                 the SAXS fits only
    benchmark-scaling-debye.svg, ...-dark.svg               the pure Debye sum only
    benchmark-accuracy.svg, benchmark-accuracy-dark.svg    deviation from the exact Debye sum, every program
    benchmark-accuracy-sas.svg, ...-dark.svg                the SAS programs, q up to 2 Å⁻¹
    benchmark-accuracy-powder.svg, ...-dark.svg             the pure Debye programs, q up to 10 Å⁻¹
    benchmark-gradient.svg, benchmark-gradient-dark.svg    time of the gradient of the Debye sum against size
    benchmark-scaling-modern.svg, ...-dark.svg              AUSAXS alone on a modern desktop, from LOCAL_RESULTS

Use each as a pair:

    ![alt](assets/figures/<name>.svg#only-light)
    ![alt](assets/figures/<name>-dark.svg#only-dark)

Points are the mean of the timed runs. A single run over the benchmark's time budget ends its series; it is drawn as
an open marker without a line. Deviations are |I/I_exact - 1| after removing each program's constant normalisation
(AUSAXS sums over ordered pairs, twice DebyeCalculator's sum). The reference is the exact Debye sum in float64, computed
here from integer squared distances (see structure_sum); no program is the judge. The SAS programs' vacuum profiles come
from docs/scripts/run_sas_accuracy.py.

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
# The accuracy of the pure Debye programs, on the silver particle; AUSAXS at its default bin width for powder diffraction.
ACCURACY_AUSAXS = "bw0.05"
ACCURACY_SIZE = "Ag_R30"
# The accuracy of the SAS programs, AUSAXS included at its SAS defaults, on a protein (see run_sas_accuracy.py).
ACCURACY_SAS_STRUCTURE = "SASDME4"
# Ag, Waasmaier & Kirfel (1995), as in include/core/form_factor/xray/XrayFormFactorTable.h; s = q/4π. DebyeCalculator and
# AES-Debye use the same parametrisation.
AG_FORM_FACTOR = ((6.073874, 17.155437, 4.173344, 0.852238, 17.988686),
                  (0.055333, 7.896512, 28.443739, 110.376106, 0.716809), 0.756603)
# The competitors run by docs/scripts/run_sas_accuracy.py: (program, label)
SAS_PROGRAMS = (("crysol", "CRYSOL"), ("pepsi-saxs", "Pepsi-SAXS"), ("foxs", "FoXS"))
# calibration molecules of run_sas_accuracy.py, by the name in the first line of their file
SAS_CALIBRATION = {"one atom": np.zeros((1, 3)), "two atoms 3 A apart": np.array([[0, 0, -1.5], [0, 0, 1.5]])}


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


def particle(folder: Path, structure: str) -> np.ndarray:
    """The coordinates of the benchmark particle, as benchmark_paper wrote them."""
    for root in (folder, LOCAL_RESULTS):
        path = root / "work" / "particles" / f"{structure}.xyz"
        if path.exists():
            return np.loadtxt(path, skiprows=2, usecols=(1, 2, 3))
    raise FileNotFoundError(f"{structure}.xyz in {folder} or {LOCAL_RESULTS}")


def structure_sum(xyz: np.ndarray, q: np.ndarray, decimals: int) -> np.ndarray:
    """The exact Debye sum of identical point scatterers, N + 2 sum_{i<j} sinc(q r_ij), in float64.

    The coordinates are written with a fixed number of decimals, so they are integers in units of 10^-decimals Å. The squared
    distances are then exact integers, and the sum runs over the distinct distances with their exact counts."""
    scale = 10**decimals
    n = np.rint(xyz*scale).astype(np.int64)
    if np.max(np.abs(n/scale - xyz)) > 1e-9:
        raise ValueError(f"coordinates are not multiples of 1e-{decimals} Å")
    key = (decimals, n.tobytes())
    if key not in _DISTANCES:
        d2 = np.concatenate([np.sum((n[i + 1:] - n[i])**2, axis=1) for i in range(len(n) - 1)] or [np.zeros(0, np.int64)])
        d2, counts = np.unique(d2, return_counts=True)
        _DISTANCES[key] = np.sqrt(d2.astype(np.float64))/scale, counts.astype(np.float64)
    r, counts = _DISTANCES[key]
    return len(n) + 2*np.array([np.dot(counts, np.sinc(x*r/np.pi)) for x in q])


_DISTANCES = {}


def ag_form_factor(q: np.ndarray) -> np.ndarray:
    a, b, c = AG_FORM_FACTOR
    s2 = (q/(4*np.pi))**2
    return sum(ai*np.exp(-bi*s2) for ai, bi in zip(a, b)) + c


def deviation(folder: Path, program: str, variant: str, structure: str, exact: dict):
    """q and |I/I_exact - 1| of a pure Debye program, with its normalisation convention (1/2, 1 or 2) removed. The exact
    intensity is f(q)^2 times the exact structure sum, as every atom is silver."""
    q, I = np.loadtxt(folder / "iq" / f"{program}_{variant}_{structure}.dat").T
    if structure not in exact or not np.array_equal(exact[structure][0], q):
        exact[structure] = q, structure_sum(particle(folder, structure), q, 6)
    ratio = I/(ag_form_factor(q)**2*exact[structure][1])
    conv = min((0.5, 1.0, 2.0), key=lambda c: abs(ratio[0]/c - 1))
    if abs(ratio[0]/conv - 1) > 0.01:
        raise ValueError(f"{program} {variant} {structure}: no known normalisation, I/I_exact = {ratio[0]:.4f}")
    return q, np.abs(ratio/conv - 1)


def sas_deviation(folder: Path, program: str, structure: str):
    """q and |I/I_exact - 1| of the vacuum profile of a SAS program, from run_sas_accuracy.py.

    The program's profile of the structure is divided by its profile of a calibration molecule, and compared with the same
    ratio of exact structure sums: as every atom is of one element, this cancels the program's form factor and
    normalisation. Each profile is taken relative to q = 0, which no expansion truncates. The exact sums use the coordinates
    exactly as the programs read them, with the three decimals of a PDB file."""
    iq = folder / "iq"
    q, I = np.loadtxt(iq / f"{program}_default_{structure}.dat").T
    qc, Ic = np.loadtxt(iq / f"{program}_default_calibration.dat").T
    with open(iq / f"{program}_default_calibration.dat") as f:
        calibration = SAS_CALIBRATION[f.readline().split(": ")[1].split(", element")[0]]
    if not np.array_equal(q, qc) or q[0] != 0:
        raise ValueError(f"{program}: the q grids differ or do not start at 0")
    xyz = np.loadtxt(folder / "work" / "sas" / f"{structure}.xyz", skiprows=2, usecols=(1, 2, 3))
    S, Sc = structure_sum(xyz, q, 3), structure_sum(calibration, q, 3)
    measured = (I/I[0])/(Ic/Ic[0])
    expected = (S/S[0])/(Sc/Sc[0])
    return q[1:], np.abs(measured/expected - 1)[1:]


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


_ACCURACY = {}
# (key, program, label, line style, running-median width in Å⁻¹) of every accuracy curve
ACCURACY_CURVES = (("aes-debye", "aes-debye", "AES-Debye", "-", 0.25),
                   ("debyecalculator", "debyecalculator", "DebyeCalculator", "--", 0.25),
                   ("ausaxs", "ausaxs", "AUSAXS", "-", 0.25),
                   ("ausaxs-sas", "ausaxs", "AUSAXS", "-", 0.05),
                   *((p, p, label, "-.", 0.05) for p, label in SAS_PROGRAMS))


def raw_deviations(folder: Path) -> dict:
    """The unsmoothed deviation curves {key: (q, deviation)}. The exact sums take minutes for the protein, so the curves are
    cached in <folder>/accuracy_curves.npz, keyed by the size and time of every input file."""
    inputs = sorted([*(folder / "iq").glob("*.dat"), *(folder / "work" / "sas").glob("*.xyz"),
                     *(p / "work" / "particles" / f"{ACCURACY_SIZE}.xyz" for p in (folder, LOCAL_RESULTS))])
    stamp = repr([(str(p), p.stat().st_size, p.stat().st_mtime_ns) for p in inputs if p.exists()])
    cache = folder / "accuracy_curves.npz"
    if cache.exists():
        with np.load(cache) as f:
            if str(f["stamp"]) == stamp:
                return {key: (f[f"{key}/q"], f[f"{key}/d"]) for key, *_ in ACCURACY_CURVES}

    exact, raw = {}, {}
    for key, variant in (("aes-debye", "res1.0"), ("debyecalculator", "float32-f0.12"), ("ausaxs", ACCURACY_AUSAXS)):
        raw[key] = deviation(folder, key, variant, ACCURACY_SIZE, exact)
    raw["ausaxs-sas"] = sas_deviation(folder, "ausaxs", ACCURACY_SAS_STRUCTURE)
    for program, _ in SAS_PROGRAMS:
        raw[program] = sas_deviation(folder, program, ACCURACY_SAS_STRUCTURE)
    np.savez(cache, stamp=stamp, **{f"{k}/{n}": v for k, (q, d) in raw.items() for n, v in (("q", q), ("d", d))})
    return raw


def accuracy_curves(folder):
    """The deviation curves, each smoothed by a running median so the spikes at the minima of I(q) don't hide the level:
    over 0.25 Å⁻¹ for the powder programs and 0.05 Å⁻¹ for the SAS programs, whose range is much shorter.
    {key: (q, deviation, program, label, line style)}"""
    if folder not in _ACCURACY:
        raw = raw_deviations(folder)
        curves = {}
        for key, program, label, ls, width in ACCURACY_CURVES:
            q, d = raw[key]
            window = max(1, int(round(width/(q[1] - q[0])))) | 1
            curves[key] = (q, rolling_median(d, window), program, label, ls)
        _ACCURACY[folder] = curves
    return _ACCURACY[folder]


def accuracy_titles(records, folder):
    """Legend titles naming the structure of each comparison."""
    atoms = {r["structure"]: r["atoms"] for r in records if r["suite"] == "accuracy"}
    with open(folder / "work" / "sas" / f"{ACCURACY_SAS_STRUCTURE}.xyz") as f:
        sas_atoms = int(f.readline())
    return {"powder": f"Silver nanoparticle, {atoms[ACCURACY_SIZE]:,} atoms",
            "sas": f"{ACCURACY_SAS_STRUCTURE}, {sas_atoms:,} atoms"}


def accuracy_axes(ax, xlim, ylim, title, xlog=False, ncol=3, loc="upper left"):
    handles, labels = ax.get_legend_handles_labels()
    ax.set_yscale("log")
    log_decades(ax.yaxis)
    if xlog:
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(ticker.FixedLocator([0.01, 0.1, 1, 10]))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(["0.01", "0.1", "1", "10"]))
        ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("q (Å⁻¹)")
    ax.set_ylabel("|I / I$_\\mathrm{exact}$ − 1|")
    ax.legend(handles[::-1], labels[::-1], loc=loc, ncol=ncol, title=title, alignment="left")


def draw_accuracy_curves(ax, c, curves, keys, labels=None):
    for key in keys:
        q, d, program, label, ls = curves[key]
        ax.plot(q, d, color=c[program], lw=1.6, ls=ls, label=(labels or {}).get(key, label))


def draw_accuracy(records, folder, c):
    """Every program on one logarithmic q axis, so both the SAS and the wide-angle range are visible. The accuracy of a
    method does not depend on the structure it is measured on, so the two comparisons share the axis; one legend each."""
    fig, ax = figure(c, size=(8, 4.4))
    fig.subplots_adjust(left=0.08, right=0.99, top=0.85, bottom=0.11)
    curves, titles = accuracy_curves(folder), accuracy_titles(records, folder)
    # the legends go above the axes, as the curves fill it from top to bottom
    for keys, title, anchor, loc in ((("crysol", "pepsi-saxs", "foxs", "ausaxs-sas"), titles["sas"], (0, 1.01), "lower left"),
                                     (("aes-debye", "debyecalculator", "ausaxs"), titles["powder"], (1, 1.01), "lower right")):
        lines = []
        for key in keys:
            q, d, program, label, ls = curves[key]
            lines += ax.plot(q, d, color=c[program], lw=1.6, ls="-." if key == "ausaxs-sas" else ls, label=label)
        ax.add_artist(ax.legend(lines[::-1], [l.get_label() for l in lines[::-1]], loc=loc, bbox_to_anchor=anchor,
                                ncol=len(lines), title=title, alignment="left", borderaxespad=0, handlelength=1.9,
                                columnspacing=0.8, handletextpad=0.5))
    accuracy_axes(ax, (0.01, 10), (1e-9, 1), None, xlog=True)
    ax.get_legend().remove()
    return fig


def draw_accuracy_sas(records, folder, c):
    """The SAS programs over their range, AUSAXS included at its SAS defaults."""
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    draw_accuracy_curves(ax, c, accuracy_curves(folder), ("crysol", "pepsi-saxs", "foxs", "ausaxs-sas"))
    accuracy_axes(ax, (0, 2), (1e-9, 1), accuracy_titles(records, folder)["sas"], ncol=4, loc="lower right")
    return fig


def draw_accuracy_powder(records, folder, c):
    """The pure Debye programs over the wide-angle range, with AUSAXS at its default bin width for powder diffraction."""
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    draw_accuracy_curves(ax, c, accuracy_curves(folder), ("aes-debye", "debyecalculator", "ausaxs"))
    accuracy_axes(ax, (0, 10), (1e-9, 1e-3), accuracy_titles(records, folder)["powder"])
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
    """AUSAXS alone on the modern desktop, with the reference machine's AUSAXS in grey. The GPU is shown at 0.25 Å only, as
    finer bin widths are currently slow on the GPU for large particles."""
    fig, ax = figure(c)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.12)
    plot(ax, c, series(local, "debye", "ausaxs", DEBYE_TIMING[0]), "ausaxs", "CPU, bin width 0.05 Å", CPU)
    plot(ax, c, series(local, "debye", "ausaxs", "forward-gpu-bw0.25"), "ausaxs", "GPU, bin width 0.25 Å", GPU)
    n = len(ax.get_legend_handles_labels()[0])
    for variant, ls, label in ((DEBYE_TIMING[0], "--", "CPU, bin width 0.05 Å"), ("forward-gpu-bw0.25", ":", "GPU, bin width 0.25 Å")):
        (x, y), _ = series(records, "debye", "ausaxs", variant)
        ax.plot(x, y, color=c["muted"], ls=ls, lw=1.0, label=label)
    handles, labels = ax.get_legend_handles_labels()
    time_axes(ax, (1.2e2, 2.5e6), (1e-4, 3e2))
    first = ax.legend(handles[:n], labels[:n], title="AUSAXS, Ryzen 5 9600X and Radeon RX 9070 XT", alignment="left",
                      loc="upper left")
    ax.add_artist(first)
    ax.legend(handles[n:], labels[n:], title="AUSAXS, i7-8700 and GTX 1070 Ti", alignment="left", loc="upper left",
              bbox_to_anchor=(0.0, 0.75))
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
                            ("benchmark-accuracy", draw_accuracy, folder), ("benchmark-accuracy-sas", draw_accuracy_sas, folder),
                            ("benchmark-accuracy-powder", draw_accuracy_powder, folder),
                            ("benchmark-gradient", draw_gradient, folder)]:
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
