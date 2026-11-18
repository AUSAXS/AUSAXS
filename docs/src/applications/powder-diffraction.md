---
title: Powder diffraction patterns of nanocrystals - AUSAXS
description: Calculate X-ray powder diffraction patterns of nanocrystals and nanoparticles from their atomic structure with the Debye scattering equation in AUSAXS.
journey: apps
hide:
  - toc
---

# Powder diffraction

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the powder diffraction pattern of a nanocrystal or nanoparticle directly from its atomic structure. A powder is a sample of randomly oriented crystallites, which is exactly what the [Debye equation](crystals.md) describes.

Since the pattern is calculated from the actual atoms of the particle, the broadening of the peaks from its finite size enters directly, as do its surface, defects and strain. For large crystallites with sharp peaks, Rietveld refinement is the better tool. The effects of the instrument, such as polarisation and instrumental broadening, are not included, and must be accounted for before comparing with a measured pattern.

The pattern is calculated as a function of $q$, which is related to the scattering angle by $q = \frac{4\pi}{\lambda}\sin\theta$.

## Example
Here we calculate the Cu Kα powder pattern of a silver nanoparticle of 201 atoms, about 2 nm across, [Ag_crystal.cif](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/Ag_crystal.cif):

=== "Command line"

    ```bash
    ausaxs debye Ag_crystal.cif --qmax 8.15 --wavelength 1.5406
    ```

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

![The powder pattern of the silver nanoparticle shows broad peaks, the strongest near 38 degrees with a shoulder near 44 degrees, followed by weaker peaks near 65, 78 and 114 degrees.](../assets/figures/powder-pattern.svg#only-light){ width="100%" }
![The powder pattern of the silver nanoparticle shows broad peaks, the strongest near 38 degrees with a shoulder near 44 degrees, followed by weaker peaks near 65, 78 and 114 degrees.](../assets/figures/powder-pattern-dark.svg#only-dark){ width="100%" }

The peaks are broad because the particle is only about 2 nm across.

Every option of `ausaxs debye` is listed in its [documentation](../documentation/cli/debye.md).
