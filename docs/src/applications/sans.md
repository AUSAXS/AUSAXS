---
title: SANS profiles of biomolecules - AUSAXS
description: "Calculate SANS profiles of proteins from atomic models with AUSAXS: neutron scattering lengths, H2O/D2O solvent contrast and hydrogen exchange."
journey: apps
hide:
  - toc
---

# Small-angle neutron scattering (SANS)

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the SANS profile of a protein or other biomolecule directly from its atomic model, and fits it to measured data. It uses the same pipeline as for [SAXS](saxs.md); only the scattering of the atoms and the solvent changes.

## What changes for neutrons
- **The atoms.** Every nucleus scatters as a point with a constant scattering length. Hydrogen scatters about half as strongly as carbon, and with the opposite sign. Hydrogens (and deuterons) therefore matter far more than in SAXS, and should be modelled explicitly whenever possible.
- **The solvent.** The solvent is usually a mixture of H₂O and D₂O, and its composition sets both the scattering of the bulk solvent and the contrast of the hydration shell. It also reaches into the molecule, as solvent-exposed hydrogens can exchange with deuterons and affect their scattering lengths.

For a deeper understanding of the physics, see e.g. [Lytje, Pedersen & Hub (2026)](https://doi.org/10.1088/1361-6633/aead1d).

## Example
Here we fit the structure of bovine serum albumin, [4f5s.pdb](https://www.sasbdb.org/media/pdb_file/SASDWB3_fit1_model1.pdb), to a SANS measurement of it in 100% D₂O, [SASDWB3.dat](https://www.sasbdb.org/data/SASDWB3/), keeping the explicit hydrogens of the structure, and plot the result:

=== "Command line"

    ```bash
    ausaxs sans 4f5s.pdb SASDWB3.dat hydrogens --keep
    ausaxs plot --folder output/sans/SASDWB3
    ```

=== "pyAUSAXS"

    ```python
    import pyausaxs as ausaxs
    import matplotlib.pyplot as plt
    ausaxs.settings.set("radiation", "neutron")
    ausaxs.settings.set("keep_hydrogens", "true")
    mol = ausaxs.create_molecule("4f5s.pdb")
    mol.hydrate()
    result = mol.fit("SASDWB3.dat")
    result.plot()
    plt.show()
    ```

![The fitted SANS profile of BSA in D₂O follows the measured data across the full q range, with residuals scattered evenly around zero.](../assets/figures/sans-fit.svg#only-light){ width="100%" }
![The fitted SANS profile of BSA in D₂O follows the measured data across the full q range, with residuals scattered evenly around zero.](../assets/figures/sans-fit-dark.svg#only-dark){ width="100%" }

For a measurement in another H₂O/D₂O mixture, pass its D₂O fraction with `solv --d2o`, e.g. `solv --d2o 0.42`; the default is 1, pure D₂O. Every option of `ausaxs sans`, including the composition of the solvent, is listed in its [documentation](../documentation/cli/sans.md).
