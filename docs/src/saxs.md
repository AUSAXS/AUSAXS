---
title: SAXS profiles of biomolecules - AUSAXS
description: How AUSAXS calculates small-angle X-ray scattering (SAXS) profiles of proteins and other biomolecules from atomic models, with an explicit hydration shell and excluded volume, and fits them to experimental SAXS data.
journey: apps
---

# Small-angle X-ray scattering (SAXS)

AUSAXS calculates the SAXS profile of a protein or other biomolecule directly from its atomic model, and fits it to measured data. The model is read from a PDB or mmCIF file.

## What goes into the profile
- **The atoms.** Every atom scatters with its own X-ray form factor. Hydrogens missing from the model are added implicitly to the atoms that carry them, in their protonation state at pH 7. Explicit hydrogens in the file can be used instead.
- **The excluded volume.** The solvent displaced by the molecule is subtracted. AUSAXS offers three published models, from a fast approximation to an explicit grid of the molecular volume.
- **The hydration shell.** A layer of water molecules is placed around the surface of the molecule. Its density is not known in advance, so its scattering power is fitted to the data.

The calculation is valid up to $q = 1$ Å⁻¹. Beyond that, the excluded volume and hydration shell models become too imprecise.

## Fitting to data
The fit adjusts the scattering power of the hydration shell. The size of the excluded volume and the solvent density can be fitted as well. The result is a fitted curve, the $\chi^2$, the fitted parameters with their uncertainties, and the model with its hydration shell written to a PDB file.

## Running it
From the command line:
```bash
ausaxs fit model.pdb data.dat
```
Or from Python, with the [pyAUSAXS](https://github.com/AUSAXS/pyAUSAXS) package:
```python
import pyausaxs as ausaxs
mol = ausaxs.create_molecule("model.pdb")
mol.hydrate()
result = mol.fit("data.dat")
```
The [SAXS fitting documentation](saxs-fitter.md) covers every option, the output files and the graphical interface. For a worked example on real data, see the [proteins tutorial](https://sastutorials.org/Proteins/Proteins.html).

The method is described in K. Lytje and J. S. Pedersen, _J. Appl. Cryst._ (2025), doi: [10.1107/S160057672500562X](https://doi.org/10.1107/S160057672500562X).
