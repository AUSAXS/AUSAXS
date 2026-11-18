---
title: Powder diffraction patterns of nanocrystals - AUSAXS
description: Calculate X-ray powder diffraction patterns of nanocrystals and nanoparticles from their atomic structure with the Debye scattering equation in AUSAXS.
journey: apps
hide:
  - toc
---

# Powder diffraction

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the powder diffraction pattern of a nanocrystal or nanoparticle directly from its atomic structure. A powder averages over every orientation of its crystallites, which is exactly what the Debye equation calculates. No assumptions about peak shapes or long-range order are needed. Larger crystallites are better handled by Rietveld refinement.

## What goes into the pattern
- **The atoms in vacuum.** Every atom scatters with its own vacuum X-ray form factor, as for [crystals and nanoparticles](crystals.md). This holds up to $q = 10$ Å⁻¹, beyond the $8.2$ Å⁻¹ reached by a Cu Kα instrument.
- **The scattering angle.** $q$ is converted to the scattering angle $2\theta$ by $q = \frac{4\pi}{\lambda}\sin\theta$.
- **What you add yourself.** The pattern is that of the particle alone. To compare it with a measured pattern, apply the corrections of your instrument, such as polarisation and peak broadening.

How to set up the calculation, and how accurate it is, is described on the [crystals and nanoparticles](crystals.md) page.

## Example
Here we calculate the Cu Kα powder pattern of a silver nanoparticle of 201 atoms, [Ag_crystal.cif](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/Ag_crystal.cif):

=== "Command line"

    ```bash
    ausaxs debye Ag_crystal.cif --qmax 8.15 --wavelength 1.5406
    ```
    The pattern is written to `output/debye/simulated/Ag_crystal/debye.dat`, as a function of both $q$ and $2\theta$.

=== "pyAUSAXS"

    ```python
    import numpy as np
    import pyausaxs as ausaxs
    import matplotlib.pyplot as plt
    ausaxs.settings.histogram(bin_width=0.02)
    wavelength = 1.5406
    q = np.linspace(0.01, 4*np.pi/wavelength, 2000)
    mol = ausaxs.create_molecule("Ag_crystal.cif")
    _, I = mol.debye_raw(q, form_factors=True)
    two_theta = np.degrees(2*np.arcsin(q*wavelength/(4*np.pi)))
    plt.plot(two_theta, I)
    plt.show()
    ```

![The powder pattern of the silver nanoparticle shows broad peaks, the strongest near 38 degrees with a shoulder near 44 degrees, followed by weaker peaks near 65, 78 and 114 degrees.](assets/figures/powder-pattern.svg#only-light){ width="100%" }
![The powder pattern of the silver nanoparticle shows broad peaks, the strongest near 38 degrees with a shoulder near 44 degrees, followed by weaker peaks near 65, 78 and 114 degrees.](assets/figures/powder-pattern-dark.svg#only-dark){ width="100%" }

The peaks are broad because the particle is only about 2 nm across. The full list of options for `ausaxs debye` is available [here](crystals.md#options).
