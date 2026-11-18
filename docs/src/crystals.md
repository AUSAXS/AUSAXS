---
title: Debye scattering of crystals and nanoparticles - AUSAXS
description: Calculate the Debye scattering of crystals, nanoparticles and arbitrary point scatterers with AUSAXS, with or without atomic form factors, from CIF, XYZ and PDB files or plain coordinate arrays.
journey: apps
hide:
  - toc
---

# Crystals and nanoparticles

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the orientationally averaged scattering of crystals, nanoparticles and any other collection of atoms or point scatterers. Unlike for [SAXS](saxs.md), there is no solvent: no excluded volume and no hydration shell.

## What goes into the profile
- **The atoms in vacuum.** Every atom scatters with its own vacuum X-ray form factor, or with its neutron scattering length.
- **Or point scatterers.** Every atom scatters as a point with its weight. This is the bare Debye equation.

The structure is read from a CIF, XYZ or PDB file, or passed as plain coordinate arrays from Python. A crystal CIF file is read as it is: no symmetry operations are applied and the unit cell is not repeated, so the file must list every atom of the particle. Non-orthogonal unit cells are not yet supported.

## Accuracy
AUSAXS sorts the pair distances into a histogram before evaluating the Debye sum, which is what makes it fast. Crystals are the hardest case for this: many pairs share exactly the same distance, so a small error in a bin position becomes a large phase error in the intensity. AUSAXS therefore places each bin at the mean of the distances it holds, instead of at its centre, and uses narrow bins of 0.02 Å. This reproduces the exact Debye sum of nanoparticles up to 8 nm across to within 0.3% of the peak heights, all the way to $q = 10$ Å⁻¹.

## Example
Here we calculate the scattering of a silver nanoparticle of 201 atoms, about 2 nm across, [Ag_crystal.cif](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/Ag_crystal.cif), up to $q = 10$ Å⁻¹:

=== "Command line"

    ```bash
    ausaxs debye Ag_crystal.cif --qmax 10
    ```
    The profile is written to `output/debye/simulated/Ag_crystal/debye.dat`.

=== "pyAUSAXS"

    ```python
    import numpy as np
    import pyausaxs as ausaxs
    import matplotlib.pyplot as plt
    ausaxs.settings.histogram(bin_width=0.02)
    q = np.linspace(0.01, 10, 2000)
    mol = ausaxs.create_molecule("Ag_crystal.cif")
    _, I = mol.debye_raw(q, form_factors=True)
    plt.semilogy(q, I)
    plt.show()
    ```

![The scattering of the silver nanoparticle falls steeply at low q, followed by broad Bragg peaks that continue out to q = 10 inverse ångström.](assets/figures/crystal-debye.svg#only-light){ width="100%" }
![The scattering of the silver nanoparticle falls steeply at low q, followed by broad Bragg peaks that continue out to q = 10 inverse ångström.](assets/figures/crystal-debye-dark.svg#only-dark){ width="100%" }

Arbitrary point scatterers are passed from Python as coordinate arrays, each with its own weight:
```python
points = ausaxs.create_molecule(x, y, z, weights)
_, I = points.debye_raw(q)
```

## Options
| Argument           | pyAUSAXS                      | Description |
|--------------------|-------------------------------|-------------|
|`--qmin`, `--qmax`  | `q`                           | The $q$-range, in Å⁻¹. |
|`--points`          | `debye_raw(q)`                | Treat every atom as a point scatterer with its weight, instead of using its form factor. |
|`--neutron`         | `settings.set("radiation", "neutron")` | Use neutron scattering lengths instead of X-ray form factors. |
|`--bin-width`       | `settings.histogram(bin_width=...)` | The width of the distance bins, in Å. The default is 0.02 Å for `ausaxs debye`; from Python, set it explicitly. |
|`--wavelength`      |                               | Also write the scattering angle $2\theta$ for this wavelength, in Å. See [powder diffraction](powder-diffraction.md). |

Given a measurement as a second file, `ausaxs debye` instead evaluates the profile at its $q$ values and fits the scale and background to it.

The [crystal scattering example](https://github.com/AUSAXS/pyAUSAXS/blob/master/examples/5_crystal_scattering.py) compares the result against the exact Debye sum and shows the effect of the bin width.
