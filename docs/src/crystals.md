---
title: Debye scattering of crystals and nanoparticles - AUSAXS
description: Calculate the Debye scattering of crystals, nanoparticles and arbitrary point scatterers with AUSAXS, with or without atomic form factors, from CIF, XYZ and PDB files or plain coordinate arrays.
---

# Crystals and nanoparticles

The Debye equation applies to any collection of scatterers, not only to biomolecules. AUSAXS calculates the orientationally averaged scattering of crystals, nanoparticles and arbitrary point scatterers, without the solvent: there is no excluded volume and no hydration shell.

## Input
- **Structure files.** CIF, XYZ and PDB files. A crystal CIF file is read as it is: AUSAXS does not apply symmetry operations or repeat the unit cell, so the file must list every atom of the particle. Non-orthogonal unit cells are not yet supported.
- **Coordinate arrays.** Any set of points, each with its own scattering weight, can be passed directly from Python.

## With or without form factors
- **Point scatterers.** Every atom scatters as a point with its weight. This is the bare Debye equation.
- **Atomic form factors.** Every atom scatters with the X-ray form factor of its element.

## Accuracy
AUSAXS sorts the pair distances into a histogram before evaluating the Debye sum, which is what makes it fast. A crystal is the hardest case for this, because many pairs share exactly the same distance, and a small error in a bin position becomes a large phase error in the intensity. AUSAXS therefore places each bin at the mean of the distances it holds, instead of at its centre. With these weighted bins and a fine bin width, the scattering of a crystal is reproduced almost exactly up to $q = 10$ Å⁻¹. From the command line, the bin width is set with `--bin-width`, and the number of bins follows from the size of the particle. From Python, it is set with `ausaxs.settings.histogram(bin_width=..., bin_count=...)`; change the two together, so that the histogram still spans the whole particle.

## Running it
From the command line:
```bash
ausaxs debye particle.cif --qmax 10
```
This writes $I(q)$ with atomic form factors to `output/debye/<structure name>/`. The options:

| Argument           | Description |
|--------------------|-------------|
|`--qmin`, `--qmax`  | The $q$-range, in Å⁻¹. |
|`--points`          | Treat every atom as a point scatterer with its weight, instead of using atomic form factors. |
|`--bin-width`       | The width of the distance bins, in Å. Lower it for large crystals. |
|`--wavelength`      | Also write the scattering angle $2\theta$ for this wavelength, in Å. See [powder diffraction](powder-diffraction.md). |

From Python, with the [pyAUSAXS](https://github.com/AUSAXS/pyAUSAXS) package, which also accepts plain coordinate arrays:
```python
import numpy as np
import pyausaxs as ausaxs

q = np.linspace(0.01, 10, 2000)
mol = ausaxs.create_molecule("crystal.cif")
_, I_points = mol.debye_raw(q)                    # point scatterers
_, I_atoms = mol.debye_raw(q, form_factors=True)  # atomic form factors

# arbitrary point scatterers with their own weights
points = ausaxs.create_molecule(x, y, z, weights)
_, I = points.debye_raw(q)
```
The [crystal scattering example](https://github.com/AUSAXS/pyAUSAXS/blob/master/examples/5_crystal_scattering.py) compares the result against the exact Debye sum and shows the effect of the bin width.
