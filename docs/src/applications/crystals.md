---
title: Debye scattering of crystals and nanoparticles - AUSAXS
description: Calculate the Debye scattering of crystals, nanoparticles and other atomic structures with AUSAXS, from CIF, XYZ and PDB files or plain coordinate arrays.
journey: apps
hide:
  - toc
---

# Crystals and nanoparticles

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the scattering of crystals, nanoparticles and any other collection of atoms with the Debye equation, using the vacuum form factor of every atom.

The Debye equation describes a sample of randomly oriented particles. A crystal therefore gives no Bragg spots, since every reflection is spread over all orientations. Instead, each set of lattice planes with spacing $d$ shows up as a peak at $q = 2\pi/d$ in the profile.

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
    ausaxs.settings.histogram(bin_width=0.05)
    q = np.linspace(0.01, 10, 2000)
    mol = ausaxs.create_molecule("Ag_crystal.cif")
    _, I = mol.debye_raw(q, form_factors=True)
    plt.semilogy(q, I)
    plt.show()
    ```

![The scattering of the silver nanoparticle falls steeply at low q, followed by broad Bragg peaks that continue out to q = 10 inverse ångström.](../assets/figures/crystal-debye.svg#only-light){ width="100%" }
![The scattering of the silver nanoparticle falls steeply at low q, followed by broad Bragg peaks that continue out to q = 10 inverse ångström.](../assets/figures/crystal-debye-dark.svg#only-dark){ width="100%" }

For a worked example in Python, see the [crystal scattering example](https://github.com/AUSAXS/pyAUSAXS/blob/master/examples/5_crystal_scattering.py). Every option of `ausaxs debye` is listed in its [documentation](../documentation/cli/debye.md).
