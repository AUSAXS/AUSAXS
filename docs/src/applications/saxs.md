---
title: SAXS profiles of biomolecules - AUSAXS
description: Calculate SAXS profiles of proteins and biomolecules from atomic models with AUSAXS, with explicit hydration shell and excluded volume, and fit SAXS data.
journey: apps
hide:
  - toc
---

# Small-angle X-ray scattering (SAXS)

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the SAXS profile of a protein or other biomolecule directly from its atomic model, and fits it to measured data.

## What goes into the profile
![The solution scattering equals the scattering of the molecule in vacuum, plus that of its hydration shell, minus that of the excluded volume.](../assets/figures/saxs-decomposition.svg#only-light){ width="100%" }
![The solution scattering equals the scattering of the molecule in vacuum, plus that of its hydration shell, minus that of the excluded volume.](../assets/figures/saxs-decomposition-dark.svg#only-dark){ width="100%" }

- **The atoms in vacuum.** Every atom scatters with its own vacuum X-ray form factor. Missing hydrogens are added implicitly to the atoms that carry them.
- **The excluded volume.** SAS only observes _variations_ in density. The average buffer density is therefore subtracted.
- **The hydration shell.** A layer of structured water molecules is placed around the surface of the molecule.

For a deeper understanding of the physics, see e.g. [Lytje, Pedersen & Hub (2026)](https://doi.org/10.1088/1361-6633/aead1d). 

## Example
Here we fit the structure of hen egg-white lysozyme, [2epe.pdb](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/2epe.pdb), to a SAXS measurement of it, [2epe.dat](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/2epe.dat), and plot the result:

=== "Command line"

    ```bash
    ausaxs saxs 2epe.pdb 2epe.dat
    ausaxs plot --folder output/saxs/2epe
    ```

=== "pyAUSAXS"

    ```python
    import pyausaxs as ausaxs
    import matplotlib.pyplot as plt
    mol = ausaxs.create_molecule("2epe.pdb")
    mol.hydrate()
    result = mol.fit("2epe.dat")
    result.plot()
    plt.show()
    ```

![The fitted SAXS profile of lysozyme follows the measured data across the full q range, with residuals scattered evenly around zero.](../assets/figures/saxs-fit.svg#only-light){ width="100%" }
![The fitted SAXS profile of lysozyme follows the measured data across the full q range, with residuals scattered evenly around zero.](../assets/figures/saxs-fit-dark.svg#only-dark){ width="100%" }

For a worked example on real data, see the [proteins tutorial](https://sastutorials.org/Proteins/Proteins.html). Every option of `ausaxs saxs` is listed in its [documentation](../documentation/cli/fit.md).
