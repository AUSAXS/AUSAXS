---
title: SAXS profiles of biomolecules - AUSAXS
description: How AUSAXS calculates small-angle X-ray scattering (SAXS) profiles of proteins and other biomolecules from atomic models, with an explicit hydration shell and excluded volume, and fits them to experimental SAXS data.
journey: apps
---

# Small-angle X-ray scattering (SAXS)

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the SAXS profile of a protein or other biomolecule directly from its atomic model, and fits it to measured data.

## What goes into the profile
![The solution scattering equals the scattering of the molecule in vacuum, plus that of its hydration shell, minus that of the excluded volume.](assets/figures/saxs-decomposition.svg#only-light){ width="100%" }
![The solution scattering equals the scattering of the molecule in vacuum, plus that of its hydration shell, minus that of the excluded volume.](assets/figures/saxs-decomposition-dark.svg#only-dark){ width="100%" }

- **The atoms in vacuum.** Every atom scatters with its own vacuum X-ray form factor. Missing hydrogens are added implicitly to the atoms that carry them.
- **The excluded volume.** SAS only observes _variations_ in density. The average buffer density is therefore subtracted.  
- **The hydration shell.** A layer of structured water molecules is placed around the surface of the molecule.

For a deeper understanding of the physics, see e.g. [Lytje, Pedersen & Hub (2026)](https://doi.org/10.1088/1361-6633/aead1d). 

## Running it
Given a measurement "data.dat", you can fit an atomic structure file, "model.pdb" like so:

=== "Command line"

    ```bash
    ausaxs fit model.pdb data.dat
    ```

=== "pyAUSAXS"

    ```python
    import pyausaxs as ausaxs
    mol = ausaxs.create_molecule("model.pdb")
    mol.hydrate()
    result = mol.fit("data.dat")
    ```

For a worked example on real data, see the [proteins tutorial](https://sastutorials.org/Proteins/Proteins.html). The full documentation
for the `saxs` program is available [here](saxs-fitter.md).