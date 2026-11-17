#!/usr/bin/env python3
"""
Figures, table and summary of benchmark_paper. Usage:

    python3 tests/benchmark/plot_paper.py output/paper_benchmark

Writes into that folder:
    summary.csv               every measurement: min, median, mean and standard deviation of the timed runs
    fig_sas.pdf               time of a full SAS fit against the number of atoms
    fig_pure_debye.pdf        (a) time of the pure Debye sum against N, with the deviation from the exact sum below;
                              (b) deviation against time at the largest accuracy size
    fig_overall.pdf           every program together; colour is the program, marker fill the kind of calculation
    table_debyecalculator.tex forward and gradient against DebyeCalculator at ~48k and ~200k atoms

Points are the mean of the timed runs, error bars one standard deviation. Deviations are the median of
|I/I_exact - 1| over q in (0, 1] and (1, 10], after removing each program's constant normalisation convention
(AUSAXS sums over all ordered pairs, which is twice DebyeCalculator's sum); the convention is checked to be one
of 1/2, 1 or 2, so a genuine scale error is not normalised away.
"""
import csv
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# One colour per program, in every figure (categorical slots 1-6 of the default validated palette); markers repeat the
# identity so it never rests on colour alone.
STYLE = {
    "ausaxs":          {"label": "AUSAXS",          "color": "#2a78d6", "marker": "o"},
    "pepsi-saxs":      {"label": "Pepsi-SAXS",      "color": "#eb6834", "marker": "s"},
    "foxs":            {"label": "FoXS",            "color": "#1baf7a", "marker": "^"},
    "crysol":          {"label": "CRYSOL",          "color": "#eda100", "marker": "D"},
    "debyecalculator": {"label": "DebyeCalculator", "color": "#e87ba4", "marker": "v"},
    "aes-debye":       {"label": "AES-Debye",       "color": "#008300", "marker": "P"},
}
INK, MUTED = "#222222", "#888888"


def load(folder):
    records = {}
    with open(folder / "results.jsonl") as f:
        for line in f:
            r = json.loads(line)
            records[r["key"]] = r  # the last record of a key wins, so a re-measurement replaces the old one
    return list(records.values())


def stats(samples):
    return {
        "min": min(samples), "median": statistics.median(samples), "mean": statistics.fmean(samples),
        "std": statistics.stdev(samples) if len(samples) > 1 else 0.0, "n": len(samples),
    }


def ok(records):
    return [r for r in records if r["status"] == "ok" and r["samples"]]


def series(records, suite, program, variant):
    """(atoms, mean, std) sorted by atoms."""
    pts = [(r["atoms"], stats(r["samples"])) for r in ok(records)
           if r["suite"] == suite and r["program"] == program and r["variant"] == variant]
    pts.sort(key=lambda p: p[0])
    return np.array([p[0] for p in pts]), np.array([p[1]["mean"] for p in pts]), np.array([p[1]["std"] for p in pts])


def style_axes(ax, xlabel, ylabel):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", color="#e6e6e6", linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def plot_series(ax, x, y, err, program, label, filled=True, dashed=False):
    if x.size == 0:
        return
    st = STYLE[program]
    ax.errorbar(x, y, yerr=err, color=st["color"], marker=st["marker"], markersize=6, linewidth=1.5,
                linestyle="--" if dashed else "-", markerfacecolor=st["color"] if filled else "white",
                markeredgewidth=1.2, capsize=2, label=label)


def write_summary(records, folder):
    with open(folder / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["suite", "program", "variant", "structure", "atoms", "status", "n", "min_s", "median_s", "mean_s", "std_s"])
        for r in sorted(records, key=lambda r: (r["suite"], r["program"], r["variant"], r["atoms"])):
            s = stats(r["samples"]) if r["samples"] else None
            w.writerow([r["suite"], r["program"], r["variant"], r["structure"], r["atoms"], r["status"],
                        s["n"] if s else 0, *(f"{s[k]:.6g}" if s else "" for k in ("min", "median", "mean", "std"))])


def fig_sas(records, folder):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    plot_series(ax, *series(records, "sas", "ausaxs-cpu", "cli"), "ausaxs", "AUSAXS (CPU)")
    plot_series(ax, *series(records, "sas", "ausaxs-gpu", "cli"), "ausaxs", "AUSAXS (GPU)", filled=False, dashed=True)
    plot_series(ax, *series(records, "sas", "pepsi-saxs", "cli"), "pepsi-saxs", "Pepsi-SAXS")
    plot_series(ax, *series(records, "sas", "foxs", "cli"), "foxs", "FoXS")
    plot_series(ax, *series(records, "sas", "crysol", "cli"), "crysol", "CRYSOL")
    style_axes(ax, "Number of atoms", "Time of a full fit (s)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(folder / "fig_sas.pdf")
    plt.close(fig)


# The normalisation convention of each program relative to the reference, i.e. I_program / I_reference at q -> 0
def deviation(iq, ref):
    q, I = iq[:, 0], iq[:, 1]
    if q.size != ref.shape[0] or np.max(np.abs(q - ref[:, 0])) > 1e-5:
        raise ValueError("q grids differ")
    ratio = I/ref[:, 1]
    conv = min((0.5, 1.0, 2.0), key=lambda c: abs(ratio[0]/c - 1))
    if abs(ratio[0]/conv - 1) > 0.01:
        raise ValueError(f"no known normalisation convention: I/I_exact = {ratio[0]:.4f} at q = {q[0]}")
    dev = np.abs(ratio/conv - 1)
    return np.median(dev[q <= 1]), np.median(dev[q > 1])


def accuracy_points(records, folder):
    """{(program, variant, structure): (atoms, mean time, dev q<=1, dev q>1)} for every profile with a reference."""
    out = {}
    iq = folder / "iq"
    for r in ok(records):
        if r["suite"] != "accuracy" or r["variant"] == "exact-float64":
            continue
        ref_file = iq / f"debyecalculator_exact-float64_{r['structure']}.dat"
        file = iq / f"{r['program']}_{r['variant']}_{r['structure']}.dat"
        if not ref_file.exists() or not file.exists():
            continue
        lo, hi = deviation(np.loadtxt(file), np.loadtxt(ref_file))
        out[(r["program"], r["variant"], r["structure"])] = (r["atoms"], stats(r["samples"])["mean"], lo, hi)
    return out


def fig_pure_debye(records, folder):
    acc = accuracy_points(records, folder)
    fig = plt.figure(figsize=(12, 5.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1], width_ratios=[1.25, 1], hspace=0.08, wspace=0.28)
    ax_t = fig.add_subplot(gs[0, 0])
    ax_a = fig.add_subplot(gs[1, 0], sharex=ax_t)
    ax_b = fig.add_subplot(gs[:, 1])

    # (a) time
    plot_series(ax_t, *series(records, "debye", "ausaxs", "forward-cpu-bw0.25"), "ausaxs", "AUSAXS, CPU")
    plot_series(ax_t, *series(records, "debye", "ausaxs", "forward-gpu-bw0.25"), "ausaxs", "AUSAXS, GPU", filled=False, dashed=True)
    for dev, filled, dashed, tag in (("cpu", True, False, "CPU"), ("cuda", False, True, "GPU")):
        plot_series(ax_t, *series(records, "debye", "debyecalculator", f"forward-{dev}"), "debyecalculator", f"DebyeCalculator, {tag}", filled, dashed)
        plot_series(ax_t, *series(records, "debye", "aes-debye", f"forward-{dev}"), "aes-debye", f"AES-Debye, {tag}", filled, dashed)
    style_axes(ax_t, "", "Time, forward (s)")
    ax_t.tick_params(labelbottom=False)
    ax_t.legend(frameon=False, fontsize=8, ncol=2)
    ax_t.set_title("(a)", loc="left", fontsize=10)

    # (a) accuracy at each program's default
    defaults = {"ausaxs": "bw0.25", "debyecalculator": "float32-f0.12", "aes-debye": "res1.0"}
    for program, variant in defaults.items():
        pts = sorted(v for (p, var, _), v in acc.items() if p == program and var == variant)
        if pts:
            st = STYLE[program]
            ax_a.plot([p[0] for p in pts], [p[3] for p in pts], color=st["color"], marker=st["marker"], linewidth=1.5, label=f"{st['label']}, $q > 1$")
            ax_a.plot([p[0] for p in pts], [p[2] for p in pts], color=st["color"], marker=st["marker"], linewidth=1.0, linestyle=":", markerfacecolor="white", label=f"{st['label']}, $q \\leq 1$")
    style_axes(ax_a, "Number of atoms", "Median deviation")
    if ax_a.lines:
        ax_a.legend(frameon=False, fontsize=7, ncol=3)

    # (b) accuracy against cost at the largest accuracy size
    sizes = sorted({v[0] for v in acc.values()})
    if sizes:
        n = sizes[-1]
        for program in ("ausaxs", "debyecalculator", "aes-debye"):
            pts = sorted(((v[1], v[3], var) for (p, var, _), v in acc.items() if p == program and v[0] == n), key=lambda p: p[0])
            if not pts:
                continue
            st = STYLE[program]
            # only AUSAXS's variants form a sweep of one parameter (the bin width); the others are unordered settings
            ax_b.plot([p[0] for p in pts], [p[1] for p in pts], color=st["color"], marker=st["marker"], markersize=7, linewidth=1.5,
                      linestyle="-" if program == "ausaxs" else "none", label=st["label"])
            for t, d, var in pts:
                ax_b.annotate(var, (t, d), textcoords="offset points", xytext=(5, 4), fontsize=7, color=MUTED)
        style_axes(ax_b, "Time, forward (s)", "Median deviation from the exact sum, $1 < q \\leq 10$")
        ax_b.set_title(f"(b) {n} atoms", loc="left", fontsize=10)
        ax_b.legend(frameon=False, fontsize=9)
    fig.savefig(folder / "fig_pure_debye.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_overall(records, folder):
    """Every program; filled markers for a fit with solvent models, open for a binned Debye sum, crossed for the exact sum."""
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for program, key in (("ausaxs", "ausaxs-cpu"), ("pepsi-saxs", "pepsi-saxs"), ("foxs", "foxs"), ("crysol", "crysol")):
        plot_series(ax, *series(records, "sas", key, "cli"), program, f"{STYLE[program]['label']}, fit")
    plot_series(ax, *series(records, "debye", "ausaxs", "forward-cpu-bw0.25"), "ausaxs", "AUSAXS, binned Debye", filled=False)
    plot_series(ax, *series(records, "debye", "debyecalculator", "forward-cpu"), "debyecalculator", "DebyeCalculator, binned Debye", filled=False)
    plot_series(ax, *series(records, "debye", "aes-debye", "forward-cpu"), "aes-debye", "AES-Debye, binned Debye", filled=False)
    x, y, e = series(records, "accuracy", "debyecalculator", "exact-float64")
    if x.size:
        st = STYLE["debyecalculator"]
        ax.plot(x, y, color=st["color"], marker="x", markersize=7, linestyle=":", label="DebyeCalculator, exact sum")
    style_axes(ax, "Number of atoms", "Time (s)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(folder / "fig_overall.pdf")
    plt.close(fig)


def table_debyecalculator(records, folder):
    """Forward and gradient, CPU, at the two particle sizes nearest 48k and 200k atoms."""
    by = defaultdict(dict)
    for r in ok(records):
        if r["suite"] == "debye":
            by[r["atoms"]][(r["program"], r["variant"])] = stats(r["samples"])["mean"]
    if not by:
        return
    rows = []
    for target in (48_000, 200_000):
        n = min(by, key=lambda a: abs(a - target))
        for mode, ours, theirs in (("forward", "forward-cpu-bw0.25", "forward-cpu"), ("gradient", "gradient-cpu", "gradient-cpu")):
            a, d = by[n].get(("ausaxs", ours)), by[n].get(("debyecalculator", theirs))
            if a and d:
                rows.append(f"{n/1000:.0f}k & {mode} & {a:.3g}\\,s & {d:.3g}\\,s & {d/a:.0f}$\\times$ \\\\")
    (folder / "table_debyecalculator.tex").write_text("\n".join(rows) + "\n")


def main():
    folder = Path(sys.argv[1] if len(sys.argv) > 1 else "output/paper_benchmark")
    records = load(folder)
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
    write_summary(records, folder)
    fig_sas(records, folder)
    fig_pure_debye(records, folder)
    fig_overall(records, folder)
    table_debyecalculator(records, folder)
    failed = [r for r in records if r["status"].startswith("failed")]
    for r in failed:
        print(f"FAILED {r['key']}: {r['status']}")
    print(f"{len(records)} measurements, {len(failed)} failed; written to {folder}")


if __name__ == "__main__":
    main()
