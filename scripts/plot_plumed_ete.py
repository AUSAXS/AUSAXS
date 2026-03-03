#!/usr/bin/env python3
"""
Visualise end-to-end distance evolution from a PLUMED MaxEnt multi-replica run.

Usage:
    python scripts/plot_plumed_ete.py <prod_dir> [--target <nm>] [--out <file.png>]

where <prod_dir> contains rep0/, rep1/, … each with a COLVAR.<i> file.

Panels:
  Top    – Per-replica ete(t) traces + ensemble average + target (time series)
  Middle – Evolving ete distribution: 2-D histogram (time × ete) for the joint
           sample of all replicas, overlaid with a rolling ensemble mean.
  Bottom – Lambda (Lagrange multiplier) evolution per replica.
"""

import argparse
import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_colvar(path: str) -> dict[str, np.ndarray]:
    """Parse a PLUMED COLVAR / LAGMULT file into a dict of named columns."""
    with open(path) as f:
        header = f.readline()
    if not header.startswith("#! FIELDS"):
        raise ValueError(f"Unexpected header in {path}: {header!r}")
    fields = header.split()[2:]  # strip '#! FIELDS'
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return {f: data[:, i] for i, f in enumerate(fields)}


def find_colvar_files(prod_dir: str) -> list[tuple[int, str]]:
    """Return sorted (replica_index, path) pairs."""
    results = []
    for path in glob.glob(os.path.join(prod_dir, "rep*", "COLVAR.*")):
        rep_dir = os.path.basename(os.path.dirname(path))
        try:
            idx = int(rep_dir.replace("rep", ""))
        except ValueError:
            continue
        results.append((idx, path))
    results.sort()
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("prod_dir",
                        help="Production directory containing rep0/, rep1/, …")
    parser.add_argument("--target", type=float, default=None,
                        help="MaxEnt target distance in nm (drawn as dashed line).")
    parser.add_argument("--out", default=None,
                        help="Save figure to this path instead of showing it.")
    parser.add_argument("--bin-dt", type=float, default=None,
                        help="Time bin width for the 2-D histogram in ps (default: auto).")

    args = parser.parse_args()

    colvar_files = find_colvar_files(args.prod_dir)
    if not colvar_files:
        sys.exit(f"No COLVAR files found under {args.prod_dir}")

    n_rep = len(colvar_files)
    colors = plt.cm.tab10(np.linspace(0, 0.9, n_rep))

    # ---- load data --------------------------------------------------------
    replicas: list[dict] = []
    for idx, path in colvar_files:
        d = load_colvar(path)
        d["_rep"] = idx
        replicas.append(d)

    time = replicas[0]["time"]  # same for all replicas
    t_max = time[-1]

    # Ensemble average is the same in every COLVAR (PLUMED broadcasts it)
    ens_ete = replicas[0]["ens_ete.ete"]

    # Lambda column name
    lambda_col = "maxent.ens_ete.ete_coupling"

    # ---- rolling ensemble mean (window ≈ 5 % of trajectory) ---------------
    win = max(1, len(time) // 20)
    rolled_mean = np.convolve(ens_ete, np.ones(win) / win, mode="same")

    # ---- figure layout ----------------------------------------------------
    fig = plt.figure(figsize=(13, 11))
    gs  = fig.add_gridspec(3, 1, hspace=0.42, height_ratios=[1.4, 1.8, 1.0])
    ax_ts  = fig.add_subplot(gs[0])   # time series
    ax_2d  = fig.add_subplot(gs[1])   # 2-D histogram
    ax_lam = fig.add_subplot(gs[2])   # lambda

    # -----------------------------------------------------------------------
    # Panel 1 – Time series
    # -----------------------------------------------------------------------
    for r, col in zip(replicas, colors):
        ax_ts.plot(r["time"], r["ete"], color=col, lw=0.8, alpha=0.7,
                   label=f"rep{r['_rep']}")

    ax_ts.plot(time, ens_ete, color="black", lw=1.8, ls="-",
               label="ensemble ⟨ete⟩", zorder=5)
    ax_ts.plot(time, rolled_mean, color="black", lw=0.9, ls="--",
               alpha=0.5, label=f"rolling mean (±{win} frames)")

    if args.target is not None:
        ax_ts.axhline(args.target, color="red", ls=":", lw=1.5,
                      label=f"target {args.target} nm")

    ax_ts.set_xlabel("Time (ps)")
    ax_ts.set_ylabel("End-to-end distance (nm)")
    ax_ts.set_title("Per-replica $r_{\\mathrm{ete}}(t)$ and ensemble average")
    ax_ts.legend(fontsize=7, ncol=2, loc="upper right")
    ax_ts.set_xlim(0, t_max)

    # -----------------------------------------------------------------------
    # Panel 2 – 2-D histogram: sampled state density vs time
    # -----------------------------------------------------------------------
    all_time = np.concatenate([r["time"] for r in replicas])
    all_ete  = np.concatenate([r["ete"]  for r in replicas])

    n_t_bins  = min(200, len(time))
    dt_bin    = args.bin_dt if args.bin_dt else t_max / n_t_bins
    ete_min   = max(0, all_ete.min() - 0.05)
    ete_max   = all_ete.max() + 0.05
    n_ete_bins = 80

    h, t_edges, e_edges = np.histogram2d(
        all_time, all_ete,
        bins=[n_t_bins, n_ete_bins],
        range=[[0, t_max], [ete_min, ete_max]]
    )
    # Normalise each time-column to a probability density
    col_sum = h.sum(axis=1, keepdims=True)
    col_sum[col_sum == 0] = 1
    h_norm = h / col_sum

    im = ax_2d.pcolormesh(
        t_edges, e_edges, h_norm.T,
        cmap="YlOrRd",
        norm=mcolors.PowerNorm(gamma=0.4),
        shading="flat",
    )
    plt.colorbar(im, ax=ax_2d, label="Normalised density (per time bin)")

    ax_2d.plot(time, ens_ete, color="black", lw=1.5, label="ensemble ⟨ete⟩")
    if args.target is not None:
        ax_2d.axhline(args.target, color="red", ls=":", lw=1.5,
                      label=f"target {args.target} nm")

    ax_2d.set_xlabel("Time (ps)")
    ax_2d.set_ylabel("End-to-end distance (nm)")
    ax_2d.set_title("Evolving sampled-state distribution (all replicas, normalised per time bin)")
    ax_2d.legend(fontsize=8, loc="upper right")
    ax_2d.set_xlim(0, t_max)
    ax_2d.set_ylim(ete_min, ete_max)

    # -----------------------------------------------------------------------
    # Panel 3 – Lambda (Lagrange multiplier) evolution
    # -----------------------------------------------------------------------
    for r, col in zip(replicas, colors):
        if lambda_col in r:
            ax_lam.plot(r["time"], r[lambda_col], color=col, lw=0.9, alpha=0.85,
                        label=f"rep{r['_rep']}")

    ax_lam.axhline(0, color="grey", ls="--", lw=0.8)
    ax_lam.set_xlabel("Time (ps)")
    ax_lam.set_ylabel("λ (kJ mol⁻¹ nm⁻¹)")
    ax_lam.set_title("MaxEnt Lagrange multiplier per replica")
    ax_lam.legend(fontsize=7, ncol=2, loc="upper left")
    ax_lam.set_xlim(0, t_max)

    # -----------------------------------------------------------------------
    # Finish
    # -----------------------------------------------------------------------
    target_str = f"  |  target = {args.target} nm" if args.target else ""
    fig.suptitle(
        f"PLUMED MaxEnt multi-replica run  ({n_rep} replicas, {t_max:.0f} ps){target_str}",
        fontsize=12, y=0.995
    )

    if args.out:
        fig.savefig(args.out, dpi=150, bbox_inches="tight")
        print(f"Saved to {args.out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
