#!/usr/bin/env python3
"""
Parse Catch2 XML benchmark output from benchmark_pr and produce plots.

Usage:
    # Run the benchmarks (from the project root):
    ./build/bin/benchmark_pr "[benchmark]" -r xml::out=tests/benchmark/results.xml

    # For a faster run with fewer samples:
    ./build/bin/benchmark_pr "[benchmark]" --benchmark-samples 10 -r xml::out=tests/benchmark/results.xml

    # Plot:
    python3 tests/benchmark/plot_benchmark.py [results.xml]
"""

import sys
import re
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

params = {
    "legend.fontsize": 16,
    "figure.figsize": (10, 7),
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "lines.markersize": 10,
    "lines.linewidth": 2,
}
plt.rcParams.update(params)

# ── helpers ──────────────────────────────────────────────────────────────────

def ns_to_ms(ns: float) -> float:
    return ns * 1e-6


def parse_atoms(section_name: str) -> int:
    """Extract atom count from e.g. 'SASDE35 (~2917 atoms)' → 2917, or 0."""
    m = re.search(r"~?(\d+)\s+atoms", section_name)
    return int(m.group(1)) if m else 0


def parse_xml(path: str) -> dict:
    """
    Returns:
        {
          "basic": {
            "Simple":   [(atoms, mean_ms, std_ms), ...],
            "Fraser":   [...],
            "Grid":     [...],
          },
          "partial": {
            # section_name → {bench_name: (mean_ms, std_ms)}
            "SASDJG5 (chain split, 2 bodies)": {"Full (baseline)": (m, s), ...},
            ...
          },
        }
    """
    tree = ET.parse(path)
    root = tree.getroot()

    basic = {}
    partial = {}  # ordered: section_name → {bench_name: (mean_ms, std_ms)}

    for tc in root.iter("TestCase"):
        name = tc.get("name", "")
        is_basic = "real molecules" in name
        is_partial = "rigidbody" in name

        for section in tc.iter("Section"):
            sec_name = section.get("name", "")

            if is_basic:
                atoms = parse_atoms(sec_name)
                if atoms == 0:
                    continue
                for br in section.findall("BenchmarkResults"):
                    bench_name = br.get("name", "")
                    mean_el = br.find("mean")
                    std_el  = br.find("standardDeviation")
                    if mean_el is None or std_el is None:
                        continue
                    mean_ms = ns_to_ms(float(mean_el.get("value")))
                    std_ms  = ns_to_ms(float(std_el.get("value")))
                    basic.setdefault(bench_name, []).append((atoms, mean_ms, std_ms))

            elif is_partial:
                bench_data = {}
                for br in section.findall("BenchmarkResults"):
                    bench_name = br.get("name", "")
                    mean_el = br.find("mean")
                    std_el  = br.find("standardDeviation")
                    if mean_el is None or std_el is None:
                        continue
                    bench_data[bench_name] = (
                        ns_to_ms(float(mean_el.get("value"))),
                        ns_to_ms(float(std_el.get("value"))),
                    )
                if bench_data:
                    partial[sec_name] = bench_data

    # sort basic series by atom count
    for k in basic:
        basic[k].sort(key=lambda x: x[0])

    return {"basic": basic, "partial": partial}


# ── plotting ─────────────────────────────────────────────────────────────────

COLOURS = {
    "Simple":          "tab:blue",
    "Fraser":          "tab:green",
    "Grid":            "tab:red",
    "Full (baseline)": "tab:orange",
    "PartialMT":       "tab:purple",
    "Partial":         "tab:brown",
}


def plot_basic(ax, data: dict, title: str):
    """Log-log scatter: atom count (x) vs time (y) for full-recalc benchmarks."""
    for label, points in data.items():
        if not points:
            continue
        pts   = np.array(points)
        colour = COLOURS.get(label, None)
        ax.errorbar(pts[:, 0], pts[:, 1], yerr=pts[:, 2],
                    fmt="o", color=colour, label=label, capsize=4)
        ax.plot(pts[:, 0], pts[:, 1], "--", color=colour)

    ax.set_xlabel("Number of atoms")
    ax.set_ylabel("Time (ms)")
    ax.set_title(title)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2)


def plot_partial(ax, data: dict, title: str):
    """Grouped bar chart: section name (x) vs time (y) for incremental benchmarks."""
    sections = list(data.keys())
    # collect all bench names across sections
    bench_names = []
    for sec_data in data.values():
        for bn in sec_data:
            if bn not in bench_names:
                bench_names.append(bn)

    x = np.arange(len(sections))
    width = 0.8 / max(len(bench_names), 1)

    for i, bench_name in enumerate(bench_names):
        means = []
        stds  = []
        for sec in sections:
            m, s = data[sec].get(bench_name, (0.0, 0.0))
            means.append(m)
            stds.append(s)
        offset = (i - (len(bench_names) - 1) / 2) * width
        colour = COLOURS.get(bench_name, None)
        ax.bar(x + offset, means, width * 0.9, yerr=stds, capsize=4,
               label=bench_name, color=colour, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(sections, rotation=15, ha="right", fontsize=13)
    ax.set_ylabel("Time (ms)")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()


def main():
    xml_path = sys.argv[1] if len(sys.argv) > 1 else "tests/benchmark/results.xml"
    if not Path(xml_path).exists():
        print(f"File not found: {xml_path}")
        print(__doc__)
        sys.exit(1)

    data = parse_xml(xml_path)

    n_plots = int(bool(data["basic"])) + int(bool(data["partial"]))
    if n_plots == 0:
        print("No benchmark data found in the XML.")
        sys.exit(1)

    fig, axes = plt.subplots(1, n_plots, figsize=(10 * n_plots, 7), squeeze=False)
    col = 0

    if data["basic"]:
        plot_basic(axes[0][col], data["basic"],
                   "Histogram managers: full recalculation")
        col += 1

    if data["partial"]:
        plot_partial(axes[0][col], data["partial"],
                     "Partial histogram managers: incremental update")

    plt.tight_layout()
    out = Path(xml_path).with_suffix(".png")
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
