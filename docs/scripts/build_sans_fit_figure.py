#!/usr/bin/env python3
"""Build the BSA fit figure shown on the SANS application page.

Fits bovine serum albumin (PDB 4f5s), keeping its explicit hydrogens, to the SEC-SANS measurement SASDWB3 in
100% D2O, exactly as in the page's example, and draws it like the SAXS fit figure (see build_fit_figure.py). Both files are
downloaded from SASBDB. The BSA drawing in the corner is figures/bsa-molecule.svg, made by build_molecule_drawing.py. Writes sans-fit.svg and sans-fit-dark.svg to docs/src/assets/figures/; use them as a pair:

    ![alt](../assets/figures/sans-fit.svg#only-light)
    ![alt](../assets/figures/sans-fit-dark.svg#only-dark)

Needs pyausaxs and matplotlib. Usage: python docs/scripts/build_sans_fit_figure.py
"""
import io
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import pyausaxs as ausaxs

from build_figures import OUTPUT
from build_fit_figure import DOCS, THEMES, draw, embed_molecule

ENTRY = "SASDWB3"
DATA_URL = f"https://www.sasbdb.org/media/intensities_files/{ENTRY}.dat"
MODEL_URL = f"https://www.sasbdb.org/media/pdb_file/{ENTRY}_fit1_model1.pdb"
MOLECULE = DOCS / "figures" / "bsa-molecule.svg"


def run_fit(folder: Path):
    data, model = folder / f"{ENTRY}.dat", folder / "4f5s.pdb"
    urllib.request.urlretrieve(DATA_URL, data)
    urllib.request.urlretrieve(MODEL_URL, model)

    ausaxs.settings.set("radiation", "neutron")
    ausaxs.settings.set("keep_hydrogens", "true")
    mol = ausaxs.create_molecule(str(model))
    mol.hydrate()
    return mol.fit(str(data))


def molecule_svg() -> tuple[str, float]:
    """The BSA drawing as an <svg> fragment with prefixed ids, and its height/width aspect."""
    text = MOLECULE.read_text(encoding="utf-8")
    _, _, width, height = map(float, re.search(r'viewBox="([^"]+)"', text).group(1).split())
    text = re.sub(r'\s(width|height)="[^"]*"', "", text, count=2)
    text = text.replace("<svg ", '<svg preserveAspectRatio="xMinYMax meet" ', 1)
    return re.sub(r'\bid="([^"]+)"', r'id="mol-\1"', text), height / width


def main() -> int:
    with tempfile.TemporaryDirectory() as folder:
        fit = run_fit(Path(folder))
    print(f"chi2r = {fit.chi2r():.3f}, parameters: {fit.fit_parameters()}")
    molecule, aspect = molecule_svg()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for scheme, colours in THEMES.items():
        fig, top = draw(fit, colours, data_label=f"{ENTRY}.dat")
        buffer = io.StringIO()
        fig.savefig(buffer, format="svg", metadata={"Date": None})
        plt.close(fig)
        out = OUTPUT / ("sans-fit.svg" if scheme == "light" else "sans-fit-dark.svg")
        out.write_text(embed_molecule(buffer.getvalue(), fig, top, molecule, aspect), encoding="utf-8")
        print(f"-> {out.relative_to(DOCS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
