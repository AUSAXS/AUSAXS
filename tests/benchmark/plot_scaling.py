#!/usr/bin/env python3
"""
Plot the results of run_scaling.py: execution time as a function of molecule size.

Usage:
    python3 tests/benchmark/plot_scaling.py [results.json ...] [-o figure.png]

With a single input this is the classic comparison figure: one dashed line per program on
log-log axes. With several inputs — e.g. the same benchmark run against two builds — every
file gets its own line style and a second panel shows the ratio to the first file, which is
the readable form of "did this change make anything faster".

    python3 tests/benchmark/plot_scaling.py before.json after.json --labels before,after

Colours are the Okabe-Ito palette, which stays separable under the common forms of colour
vision deficiency. Identity is never carried by colour alone: every series also has its own
marker, and the AUSAXS series are drawn with filled markers against the open markers of the
other programs.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

params = {
    "legend.fontsize": 15,
    "figure.figsize": (10, 8),
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "lines.markersize": 9,
    "lines.linewidth": 2,
}
plt.rcParams.update(params)

# Fixed assignment, never cycled: a series keeps its colour no matter which subset was run.
# 'ausaxs' marks the series drawn with a filled marker.
STYLE = {
    "pepsi":            {"label": "Pepsi-SAXS",                          "color": "#0072B2", "marker": "o"},
    "foxs":             {"label": "FoXS",                                "color": "#E69F00", "marker": "s"},
    "crysol":           {"label": "CRYSOL",                              "color": "#56B4E9", "marker": "^"},
    "ausaxs_simple":    {"label": "AUSAXS$_{simple}$",                   "color": "#D55E00", "marker": "o", "ausaxs": True},
    "ausaxs_simple_gpu": {"label": "AUSAXS$_{simple}$ (GPU)",            "color": "#2E2E2E", "marker": "P", "ausaxs": True},
    "ausaxs_simple_st": {"label": "AUSAXS$_{simple}$ (1 thread)",        "color": "#6A3D9A", "marker": "v", "ausaxs": True},
    "ausaxs_fraser":    {"label": "AUSAXS$_{gaussian\\ spheres}$",         "color": "#009E73", "marker": "D", "ausaxs": True},
    "ausaxs_grid":      {"label": "AUSAXS$_{grid}$",               "color": "#CC79A7", "marker": "s", "ausaxs": True},
}
ORDER = list(STYLE)
LINESTYLES = ["--", "-", ":", "-."]


def load(path: Path) -> dict:
    """{method: (atoms, mean_ms, std_ms) sorted by atoms}"""
    records = json.loads(path.read_text())["results"]
    series = {}
    for r in records:
        series.setdefault(r["method"], []).append((r["atoms"], r["mean"] * 1e3, r["stddev"] * 1e3))
    return {m: np.array(sorted(v)) for m, v in series.items()}


def plot_absolute(ax, datasets: list, labels: list):
    for i, (data, file_label) in enumerate(zip(datasets, labels)):
        ls = LINESTYLES[i % len(LINESTYLES)]
        for method in ORDER:
            if method not in data:
                continue
            style = STYLE[method]
            pts = data[method]
            label = style["label"] + (f" [{file_label}]" if len(datasets) > 1 else "")
            ax.errorbar(
                pts[:, 0], pts[:, 1], yerr=pts[:, 2],
                fmt=style["marker"], color=style["color"], label=label, capsize=3,
                markerfacecolor=style["color"] if style.get("ausaxs") else "none",
            )
            ax.plot(pts[:, 0], pts[:, 1], ls, color=style["color"])

    ax.set_ylabel("Execution time (ms)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(ncol=2)


def plot_ratio(ax, datasets: list, labels: list):
    """Time relative to the first dataset, on the atom counts the two have in common."""
    reference = datasets[0]
    for i, (data, file_label) in enumerate(zip(datasets[1:], labels[1:]), start=1):
        ls = LINESTYLES[i % len(LINESTYLES)]
        for method in ORDER:
            if method not in data or method not in reference:
                continue
            style = STYLE[method]
            ref = {int(a): t for a, t, _ in reference[method]}
            shared = np.array([(a, t / ref[int(a)]) for a, t, _ in data[method] if int(a) in ref])
            if shared.size == 0:
                continue
            ax.plot(shared[:, 0], shared[:, 1], ls, marker=style["marker"], color=style["color"],
                    markerfacecolor=style["color"] if style.get("ausaxs") else "none",
                    label=f"{style['label']} [{file_label}]")

    ax.axhline(1.0, color="0.4", lw=1)
    ax.set_ylabel(f"Time / {labels[0]}")
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(ncol=2, fontsize=12)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("results", nargs="*", default=["output/scaling_benchmark/results.json"],
                   help="one or more JSON files written by run_scaling.py")
    p.add_argument("-o", "--out", default="", help="output image (default: alongside the first input)")
    p.add_argument("--labels", default="", help="comma-separated names for the inputs")
    p.add_argument("--methods", default="", help="comma-separated subset of series to draw")
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args()

    paths = [Path(r) for r in args.results]
    for path in paths:
        if not path.exists():
            raise SystemExit(f"file not found: {path}. Run tests/benchmark/run_scaling.py first.")

    datasets = [load(path) for path in paths]
    if args.methods:
        wanted = set(args.methods.split(","))
        unknown = wanted - set(STYLE)
        if unknown:
            raise SystemExit(f"unknown series: {', '.join(sorted(unknown))}")
        datasets = [{m: v for m, v in d.items() if m in wanted} for d in datasets]
    labels = args.labels.split(",") if args.labels else [path.stem for path in paths]
    if len(labels) != len(paths):
        raise SystemExit("--labels must give one name per input file")

    if len(datasets) > 1:
        fig, (ax, ax_ratio) = plt.subplots(
            2, 1, figsize=(10, 11), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
        plot_absolute(ax, datasets, labels)
        plot_ratio(ax_ratio, datasets, labels)
        ax_ratio.set_xlabel("Number of atoms")
    else:
        fig, ax = plt.subplots()
        plot_absolute(ax, datasets, labels)
        ax.set_xlabel("Number of atoms")

    fig.tight_layout()
    out = Path(args.out) if args.out else paths[0].with_suffix(".png")
    fig.savefig(out, dpi=args.dpi)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
