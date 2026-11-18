---
title: SANS profiles of biomolecules - AUSAXS
description: How AUSAXS calculates small-angle neutron scattering (SANS) profiles of proteins from atomic models, with neutron scattering lengths, H2O/D2O solvent contrast and hydrogen-deuterium exchange.
journey: apps
hide:
  - toc
---

# Small-angle neutron scattering (SANS)

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the SANS profile of a protein or other biomolecule directly from its atomic model, and fits it to measured data. It uses the same pipeline as for [SAXS](saxs.md); only the scattering of the atoms and the solvent changes.

## What changes for neutrons
- **The atoms.** Every nucleus scatters as a point with a constant scattering length. Hydrogen scatters about as strongly as carbon, but with the opposite sign. Hydrogens therefore matter far more than in SAXS, and should be modeled explicitly whenever possible.
- **The solvent.** The solvent is usually a mixture of H₂O and D₂O, and its composition sets both the scattering of the bulk solvent and the contrast of the hydration shell. It also reaches into the molecule, as solvent-exposed hydrogens can exchange with deuterons and affect their scattering lengths.

For a deeper understanding of the physics, see e.g. [Lytje, Pedersen & Hub (2026)](https://doi.org/10.1088/1361-6633/aead1d).

## Example
Given a SANS measurement in 42% D₂O, "data.dat", you can fit an atomic structure file, "model.pdb", like so:

=== "Command line"

    ```bash
    ausaxs sans model.pdb data.dat solv --d2o 0.42
    ausaxs plot --folder output/sans_fitter/data
    ```

=== "pyAUSAXS"

    ```python
    import pyausaxs as ausaxs
    import matplotlib.pyplot as plt
    ausaxs.settings.set("radiation", "neutron")
    ausaxs.settings.set("neutron_solvent_density", "0.2348") # 42% D₂O
    mol = ausaxs.create_molecule("model.pdb")
    mol.hydrate()
    result = mol.fit("data.dat")
    result.plot()
    plt.show()
    ```

Every option of `ausaxs sans`, including the composition of the solvent, is listed in its [documentation](../documentation/cli/sans.md).
