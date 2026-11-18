---
title: Debye scattering of crystals and nanoparticles - AUSAXS
description: Calculate the Debye scattering of crystals, nanoparticles and arbitrary point scatterers with AUSAXS, with or without atomic form factors, from CIF, XYZ and PDB files or plain coordinate arrays.
journey: apps
hide:
  - toc
---

# Crystals and nanoparticles

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the orientationally averaged scattering of crystals, nanoparticles and any other collection of atoms or point scatterers.

## What goes into the profile
- **The atoms in vacuum.** Every atom scatters with its own vacuum X-ray form factor, or with its neutron scattering length.
- **Or point scatterers.** Every scatterer is a point with a constant weight. This is the bare Debye equation, and applies to any system that can be described as a set of points.
- **The orientational average.** The Debye equation sums over all pairs of scatterers, and only their distances enter. It therefore describes a particle observed in every orientation at once, and makes no assumption of long-range order. The Bragg peaks of a crystal appear naturally, broadened by its finite size.

## Example
Here we calculate the scattering of a silver nanoparticle of 201 atoms, about 2 nm across, [Ag_crystal.cif](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/Ag_crystal.cif), up to $q = 10$ Å⁻¹:

=== "Command line"

    ```bash
    ausaxs debye Ag_crystal.cif --qmax 10
    ```

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

![The scattering of the silver nanoparticle falls steeply at low q, followed by broad Bragg peaks that continue out to q = 10 inverse ångström.](../assets/figures/crystal-debye.svg#only-light){ width="100%" }
![The scattering of the silver nanoparticle falls steeply at low q, followed by broad Bragg peaks that continue out to q = 10 inverse ångström.](../assets/figures/crystal-debye-dark.svg#only-dark){ width="100%" }

For a worked example in Python, see the [crystal scattering example](https://github.com/AUSAXS/pyAUSAXS/blob/master/examples/5_crystal_scattering.py).
