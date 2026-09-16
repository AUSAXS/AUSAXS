#!/usr/bin/env python3
"""
End-to-end scaling benchmark: wall-clock time of a complete SAXS fit as a function of molecule size.

This is the whole-program counterpart to benchmark_pr.cpp. Where that one measures individual
library calls, this one measures what a user actually waits for: process start, file parsing,
hydration, histogram, fit and output. That is the only granularity in which AUSAXS can be
compared against other SAXS programs, so Pepsi-SAXS, FoXS and CRYSOL are included whenever they
are found on the system.

Usage (from the repository root, so the data/ paths resolve):
    cmake -B build-release -S . -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=OFF
    cmake --build build-release --target ausaxs -j8
    python3 tests/benchmark/run_scaling.py
    python3 tests/benchmark/plot_scaling.py

Useful options:
    --methods ausaxs_simple,ausaxs_grid     only run a subset
    --molecules 1ubq,SASDJG5                only run a subset
    --max-atoms 20000                       skip everything larger
    --budget 20                             seconds to spend per (molecule, method)
    --ausaxs path/to/ausaxs --out other.json    benchmark a different build (A/B comparisons)

Requires hyperfine (https://github.com/sharkdp/hyperfine).

Notes on methodology:
  * Every structure is stripped to heavy ATOM records only and copied, together with its
    measurement file, into a private working directory. Stripping puts the programs on equal
    footing (they disagree on what to do with explicit hydrogens), and the copy keeps the
    per-dataset settings.txt files living in data/ from being auto-discovered.
  * The number of runs is chosen from a single probe run against a time budget, so small systems
    are measured many times and large ones only a couple. The probe also warms the file cache.
  * The AUSAXS commands mirror the ones used for the original comparison plot: the simple model
    runs with its defaults, while the Fraser and grid models fit the excluded volume, which is
    how they are meant to be used.
"""

import argparse
import json
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ── molecules ────────────────────────────────────────────────────────────────
# A roughly log-spaced ladder through data/. Atom counts in the comments are heavy atoms after
# stripping, which is what ends up on the x-axis; they are recounted at run time.
MOLECULES = [
    ("1ubq",     "data/1ubq/1ubq.pdb",                         "data/1ubq/1ubq.dat"),                          #   602
    ("6lyz",     "data/6lyz/6lyz.pdb",                         "data/6lyz/6lyz.dat"),                          #  1001
    ("SASDPS4",  "data/consensus/old/SASDPS4/SASDPS4.pdb",     "data/consensus/old/SASDPS4/SASDPS4.dat"),      #  1481
    ("SASDAG2",  "data/SASDAG2/SASDAG2.pdb",                   "data/SASDAG2/SASDAG2.dat"),                    #  2002
    ("SASDE35",  "data/SASDE35/SASDE35.pdb",                   "data/SASDE35/SASDE35.dat"),                    #  2917
    ("SASDQ59",  "data/rigidbody/SASDQ59/SASDQ59.pdb",         "data/rigidbody/SASDQ59/SASDQ59.dat"),          #  3540
    ("SASDJG5",  "data/SASDJG5/SASDJG5.pdb",                   "data/SASDJG5/SASDJG5.dat"),                    #  4734
    ("SASDME4",  "data/SASDME4/SASDME4.pdb",                   "data/SASDME4/SASDME4.dat"),                    #  5380
    ("SASDPB9",  "data/rigidbody/SASDPB9/SASDPB9.pdb",         "data/rigidbody/SASDPB9/SASDPB9.dat"),          #  6210
    ("SHOC2",    "data/SHOC2/SHOC2.pdb",                       "data/SHOC2/SHOC2.dat"),                        #  7694
    ("urateox",  "data/rigidbody/urateox/urateox.pdb",         "data/rigidbody/urateox/urateox.dat"),          #  9436
    ("SASDDD3",  "data/SASDDD3/SASDDD3.pdb",                   "data/SASDDD3/SASDDD3.dat"),                    # 10775
    ("SASDPR4",  "data/consensus/old/SASDPR4/SASDPR4.pdb",     "data/consensus/old/SASDPR4/SASDPR4.dat"),      # 12332
    ("SASDA92",  "data/SASDA92/SASDA92.pdb",                   "data/SASDA92/SASDA92.dat"),                    # 16068
    ("SASDBC3",  "data/symmetry/SASDBC3/SASDBC3.pdb",          "data/symmetry/SASDBC3/SASDBC3.dat"),           # 19872
    ("SASDA45",  "data/SASDA45/SASDA45.pdb",                   "data/SASDA45/SASDA45.dat"),                    # 25761
    ("SASDGL2",  "data/symmetry/SASDGL2/SASDGL2.pdb",          "data/symmetry/SASDGL2/SASDGL2.dat"),           # 30564
    ("A2M_ma",   "data/A2M_ma/A2M_ma.pdb",                     "data/A2M_ma/A2M_ma.dat"),                      # 43564
]

# ── external programs ────────────────────────────────────────────────────────
# Only the ones that are actually present are benchmarked; everything else is silently skipped.
EXTERNAL_CANDIDATES = {
    "pepsi":  ["Pepsi-SAXS", str(Path.home() / "tools/Pepsi-SAXS/Pepsi-SAXS")],
    "foxs":   ["foxs",       str(Path.home() / "tools/imp/build/bin/foxs")],
    "crysol": ["crysol",     str(Path.home() / "tools/ATSAS/bin/crysol")],
}

AUSAXS_COMMON = ["--output", "out/", "--allow-unknown-atoms", "--allow-unknown-residues", "--offline"]

def ausaxs_methods(ausaxs: str) -> dict:
    return {
        "ausaxs_simple":    [ausaxs, "fit", "{pdb}", "{dat}"] + AUSAXS_COMMON,
        "ausaxs_simple_st": [ausaxs, "fit", "{pdb}", "{dat}"] + AUSAXS_COMMON + ["-t", "1"],
        "ausaxs_fraser":    [ausaxs, "fit", "{pdb}", "{dat}"] + AUSAXS_COMMON + ["exv", "--model", "fraser", "--fit"],
        "ausaxs_grid":      [ausaxs, "fit", "{pdb}", "{dat}"] + AUSAXS_COMMON + ["exv", "--model", "grid", "--fit"],
        # GPU acceleration currently covers the simple model only. It needs the backend library to
        # have been built (cmake --build <build> --target ausaxs_gpu_sycl), and falls back to the
        # CPU kernel without failing when there is none, which GPU_FALLBACK_MARKER catches below.
        "ausaxs_simple_gpu": [ausaxs, "fit", "{pdb}", "{dat}"] + AUSAXS_COMMON + ["--gpu"],
    }

GPU_METHODS = {"ausaxs_simple_gpu"}
GPU_FALLBACK_MARKER = "no usable GPU backend"

def external_methods() -> dict:
    """The external programs that are both installed and actually able to start. A binary left
    behind by an older distribution typically still exists but dies in the dynamic loader, which
    is indistinguishable from "not installed" as far as this benchmark is concerned."""
    found = {}
    for name, candidates in EXTERNAL_CANDIDATES.items():
        exe = next((c for c in candidates if shutil.which(c)), None)
        if exe is None:
            continue
        if not _can_launch(exe):
            print(f"[skip] {name}: {exe} exists but does not run (missing shared libraries?)")
            continue
        if name == "pepsi":
            found[name] = [exe, "{pdb}", "{dat}", "-o", "out.fit"]
        elif name == "foxs":
            found[name] = [exe, "{pdb}", "{dat}"]
        elif name == "crysol":
            found[name] = [exe, "{dat}", "{pdb}", "--constant", "--implicit-hydrogen", "1"]
    return found

def _can_launch(exe: str) -> bool:
    """Ask the dynamic loader whether the binary could start, rather than starting it: several of
    these programs drop into an interactive prompt when run without arguments."""
    try:
        r = subprocess.run(["ldd", exe], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return True  # no way to tell; let the probe run decide
    return "not found" not in r.stdout

# ── structure preparation ────────────────────────────────────────────────────
_WATER_RESIDUES = {"HOH", "WAT", "SOL", "DOD", "TIP", "TIP3"}

def _is_hydrogen(line: str) -> bool:
    """True for hydrogen/deuterium. The element column is empty in several of the data/ files, so
    fall back on the atom name, whose leading digits are the (optional) branch index."""
    element = line[76:78].strip()
    if element:
        return element in ("H", "D")
    name = line[12:16].strip().lstrip("0123456789")
    return name[:1] in ("H", "D")

def strip_structure(src: Path, dst: Path) -> int:
    """Write src to dst keeping only heavy, non-solvent ATOM records. Returns the atom count."""
    kept = []
    for line in src.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[17:20].strip() in _WATER_RESIDUES or line[12:16].strip() == "OW":
            continue
        if _is_hydrogen(line):
            continue
        kept.append(line)
    dst.write_text("\n".join(kept) + "\nEND\n")
    return len(kept)

# ── measurement ──────────────────────────────────────────────────────────────
def probe(cmd: list, cwd: Path) -> tuple:
    """One untimed-for-the-record run: picks the sample count below and warms the file cache.
    Returns its duration and its output, the latter only so the GPU fallback can be spotted."""
    t0 = time.perf_counter()
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, errors="replace")
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        raise RuntimeError(f"command failed with exit code {r.returncode}: {' '.join(cmd)}")
    return dt, r.stdout + r.stderr

def measure(cmd: list, cwd: Path, name: str, runs: int, out_json: Path) -> dict:
    """hyperfine --shell=none, so the reported times contain no shell spawn. Everything after the
    `--` is one benchmarked command, so the argv list has to be joined back into a single string."""
    hf = [
        "hyperfine", "--shell=none", "--style", "none",
        "--warmup", "1" if runs >= 5 else "0",
        "--runs", str(runs),
        "--command-name", name,
        "--export-json", str(out_json),
        "--", shlex.join(cmd),
    ]
    subprocess.run(hf, cwd=cwd, check=True, stdout=subprocess.DEVNULL)
    return json.loads(out_json.read_text())["results"][0]

# ── driver ───────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ausaxs", default="build-release/bin/ausaxs", help="AUSAXS binary to benchmark")
    p.add_argument("--workdir", default="output/scaling_benchmark", help="scratch directory for the prepared inputs")
    p.add_argument("--out", default="output/scaling_benchmark/results.json", help="where to write the results")
    p.add_argument("--budget", type=float, default=20.0, help="seconds to spend per (molecule, method)")
    p.add_argument("--min-runs", type=int, default=2)
    p.add_argument("--max-runs", type=int, default=10)
    p.add_argument("--max-atoms", type=int, default=0, help="skip molecules above this size (0 = no limit)")
    p.add_argument("--methods", default="", help="comma-separated subset of method keys")
    p.add_argument("--molecules", default="", help="comma-separated subset of molecule labels")
    p.add_argument("--no-external", action="store_true", help="benchmark AUSAXS only")
    p.add_argument("--force", action="store_true", help="re-measure entries already present in --out")
    args = p.parse_args()

    ausaxs = Path(args.ausaxs)
    if not ausaxs.exists():
        sys.exit(f"AUSAXS binary not found: {ausaxs}. Build it first (see the module docstring).")
    if shutil.which("hyperfine") is None:
        sys.exit("hyperfine not found on PATH; it is required to run this benchmark.")

    methods = ausaxs_methods(str(ausaxs.resolve()))
    if not args.no_external:
        methods.update(external_methods())
    if args.methods:
        wanted = args.methods.split(",")
        missing = [m for m in wanted if m not in methods]
        if missing:
            sys.exit(f"unknown/unavailable methods: {', '.join(missing)} (have: {', '.join(methods)})")
        methods = {k: v for k, v in methods.items() if k in wanted}

    molecules = MOLECULES
    if args.molecules:
        wanted = args.molecules.split(",")
        molecules = [m for m in molecules if m[0] in wanted]
        missing = set(wanted) - {m[0] for m in molecules}
        if missing:
            sys.exit(f"unknown molecules: {', '.join(sorted(missing))}")

    # absolute: the benchmarked commands run with their own directory as the working directory
    workdir = Path(args.workdir).resolve()
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # results are keyed by (molecule, method) so an interrupted run can simply be restarted
    results = {}
    if out_path.exists() and not args.force:
        results = {(r["molecule"], r["method"]): r for r in json.loads(out_path.read_text())["results"]}
        print(f"resuming from {out_path} with {len(results)} existing measurements")

    def flush():
        records = sorted(results.values(), key=lambda r: (r["atoms"], r["method"]))
        out_path.write_text(json.dumps({
            "ausaxs": str(ausaxs.resolve()),
            "host": {"cpu": _cpu_name(), "threads": _thread_count()},
            "results": records,
        }, indent=1))

    print(f"methods: {', '.join(methods)}")
    for label, pdb_src, dat_src in molecules:
        pdb_src, dat_src = Path(pdb_src), Path(dat_src)
        if not pdb_src.exists() or not dat_src.exists():
            print(f"[skip] {label}: missing input files")
            continue

        base = workdir / label
        base.mkdir(parents=True, exist_ok=True)
        stripped = base / f"{label}_stripped.pdb"
        atoms = strip_structure(pdb_src, stripped)
        if args.max_atoms and atoms > args.max_atoms:
            print(f"[skip] {label}: {atoms} atoms exceeds --max-atoms")
            continue
        print(f"\n=== {label}: {atoms} atoms ===")

        for method, template in methods.items():
            if (label, method) in results and not args.force:
                print(f"  {method:<19} cached")
                continue

            # every program gets its own directory holding a private copy of both inputs, so the
            # output files they scatter around cannot interfere with each other
            cwd = base / method
            cwd.mkdir(exist_ok=True)
            shutil.copy(stripped, cwd / stripped.name)
            shutil.copy(dat_src, cwd / dat_src.name)
            cmd = [c.format(pdb=stripped.name, dat=dat_src.name) for c in template]

            try:
                probe_t, probe_out = probe(cmd, cwd)
            except RuntimeError as e:
                print(f"  {method:<19} FAILED: {e}")
                continue
            if method in GPU_METHODS and GPU_FALLBACK_MARKER in probe_out:
                # timing a silent CPU fallback under a GPU label would be worse than no data point
                print(f"  {method:<19} skipped: fell back to the CPU kernel, no usable GPU backend")
                continue
            runs = max(args.min_runs, min(args.max_runs, int(args.budget / max(probe_t, 1e-3))))
            print(f"  {method:<19} probe {probe_t*1e3:8.1f} ms → {runs} runs", end="", flush=True)

            r = measure(cmd, cwd, method, runs, cwd / "hyperfine.json")
            results[(label, method)] = {
                "molecule": label, "method": method, "atoms": atoms, "runs": runs,
                "mean": r["mean"], "stddev": r["stddev"] or 0.0, "min": r["min"], "max": r["max"],
                "command": " ".join(cmd),
            }
            print(f"  →  {r['mean']*1e3:9.1f} ± {(r['stddev'] or 0.0)*1e3:.1f} ms")
            flush()

    flush()
    print(f"\nwrote {len(results)} measurements to {out_path}")

def _cpu_name() -> str:
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name"):
            return line.split(":", 1)[1].strip()
    return "unknown"

def _thread_count() -> int:
    import os
    return os.cpu_count() or 0

if __name__ == "__main__":
    main()
