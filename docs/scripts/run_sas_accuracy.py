#!/usr/bin/env python3
"""Run CRYSOL, Pepsi-SAXS, FoXS and AUSAXS on a protein for the SAS accuracy figure of build_benchmark_figures.py.

The SAS programs expand the scattering amplitude or bin the distances in their own way, and none can be run as a plain
Debye sum. Their accuracy is therefore measured on the vacuum (atomic) profile, which each writes beside the total profile:
CRYSOL's .int column 3, Pepsi-SAXS's Iat and FoXS's partial-profile column 2. AUSAXS runs in-process with the vacuum model.

The structure is a protein, since a SAS program is meant for a disordered molecule in solution; a crystal's few discrete
distances are a worst case for binning that a protein never presents. To measure the error against an exact sum, every
heavy atom of the protein becomes the same element, written as HETATM records so that no program adds implicit hydrogens:
silver, except for Pepsi-SAXS, which has no silver and gets zinc. Only the geometry matters here.

Each program runs with its own defaults, except that q_max is raised to its limit: 2 Å⁻¹ for CRYSOL and FoXS and 1 Å⁻¹ for
Pepsi-SAXS. AUSAXS is evaluated at the same q as CRYSOL.

Beside the protein, each program also gets a calibration molecule with an exactly known structure sum: a single atom, or
for CRYSOL, which fails on one atom, two atoms 3 Å apart. Divided by its own calibration, a program's profile of the protein
carries neither its form factor nor its normalisation, as every atom is of one element.

Writes, into <output>:

    work/sas/<label>.xyz                        the heavy-atom coordinates every program read, with three decimals
    iq/<program>_default_<label>.dat            q and the vacuum intensity of the protein
    iq/<program>_default_calibration.dat        the same for the calibration molecule; the first line says which

Usage: python docs/scripts/run_sas_accuracy.py <structure.pdb> [output/paper_benchmark]
with AUSAXS_BENCH_CRYSOL, AUSAXS_BENCH_PEPSI, AUSAXS_BENCH_FOXS and their _LIBS, as for benchmark_paper, and pyausaxs.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

DOCS = Path(__file__).resolve().parent.parent

# The calibration molecules, by name: coordinates in Å.
CALIBRATION = {"one atom": np.zeros((1, 3)), "two atoms 3 A apart": np.array([[0, 0, -1.5], [0, 0, 1.5]])}


def read_heavy_atoms(path: Path) -> np.ndarray:
    """The heavy ATOM records of the first model, as benchmark_paper gives them to every SAS program."""
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[76:78].strip() not in ("H", "D"):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def write_pdb(path: Path, element: str, xyz: np.ndarray):
    el = element.upper()
    with open(path, "w") as f:
        for i, (x, y, z) in enumerate(xyz, start=1):
            f.write(f"HETATM{i % 100000:5d} {el:<4s} {el:>3s} A{i % 10000:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
        f.write("END\n")


def run(program: str, args: list[str], cwd: Path):
    exe = os.environ.get(f"AUSAXS_BENCH_{program}")
    if not exe:
        raise RuntimeError(f"AUSAXS_BENCH_{program} is not set")
    env = dict(os.environ)
    if libs := os.environ.get(f"AUSAXS_BENCH_{program}_LIBS"):
        env["LD_LIBRARY_PATH"] = libs + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    subprocess.run([exe, *args], cwd=cwd, env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def crysol(pdb: Path) -> np.ndarray:
    run("CRYSOL", [pdb.name, "--smax=2", "--ns=1001", "-p", "out"], pdb.parent)
    d = np.loadtxt(pdb.parent / "out.int", skiprows=1)
    return d[:, [0, 2]]


def pepsi(pdb: Path) -> np.ndarray:
    run("PEPSI", [pdb.name, "-o", "out.dat", "-ms", "1.0", "-ns", "1001"], pdb.parent)
    d = np.loadtxt(pdb.parent / "out.dat", comments="#")
    return d[:, [0, 2]]


def foxs(pdb: Path) -> np.ndarray:
    run("FOXS", [pdb.name, "-p", "-q", "2", "-s", "1000"], pdb.parent)
    d = np.loadtxt(pdb.parent / f"{pdb.name}.dat", comments="#")
    return d[:, [0, 1]]


def ausaxs(pdb: Path) -> np.ndarray:
    import pyausaxs
    from pyausaxs import ExvModel

    pyausaxs.settings.exv(ExvModel.vacuum)
    q = np.linspace(0, 2, 1001)
    _, I = pyausaxs.create_molecule(str(pdb)).debye(q)
    return np.column_stack([q, I])


PROGRAMS = {
    # name: (run, element, calibration molecule)
    "crysol": (crysol, "Ag", "two atoms 3 A apart"),
    "pepsi-saxs": (pepsi, "Zn", "one atom"),
    "foxs": (foxs, "Ag", "one atom"),
    "ausaxs": (ausaxs, "Ag", "one atom"),
}


def main():
    structure = Path(sys.argv[1])
    folder = Path(sys.argv[2]) if len(sys.argv) > 2 else DOCS.parent / "output" / "paper_benchmark"
    label = structure.stem
    xyz = read_heavy_atoms(structure)
    (folder / "iq").mkdir(parents=True, exist_ok=True)
    (folder / "work" / "sas").mkdir(parents=True, exist_ok=True)
    with open(folder / "work" / "sas" / f"{label}.xyz", "w") as f:
        f.write(f"{len(xyz)}\nheavy atoms of {structure.name}, as read by the SAS programs\n")
        for x, y, z in xyz:
            f.write(f"X {x:.3f} {y:.3f} {z:.3f}\n")

    for name, (calc, element, calibration) in PROGRAMS.items():
        for kind, coords, out in (("structure", xyz, f"{name}_default_{label}.dat"),
                                  ("calibration", CALIBRATION[calibration], f"{name}_default_calibration.dat")):
            with tempfile.TemporaryDirectory() as tmp:
                pdb = Path(tmp) / "structure.pdb"
                write_pdb(pdb, element, coords)
                profile = calc(pdb)
            header = f"{name}, {kind}: {calibration if kind == 'calibration' else label}, element {element}\nq I_vacuum"
            np.savetxt(folder / "iq" / out, profile, header=header)
            print(f"-> {folder / 'iq' / out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
