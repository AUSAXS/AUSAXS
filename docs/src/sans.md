---
title: SANS profiles of biomolecules - AUSAXS
description: How AUSAXS calculates small-angle neutron scattering (SANS) profiles of proteins from atomic models, with neutron scattering lengths, H2O/D2O solvent contrast and hydrogen-deuterium exchange.
journey: apps
---

# Small-angle neutron scattering (SANS)

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> calculates the SANS profile of a protein or other biomolecule directly from its atomic model, and fits it to measured data. It uses the same pipeline as for [SAXS](saxs.md); only the scattering of the atoms and the solvent changes.

## What changes for neutrons
- **The atoms.** Every nucleus scatters as a point with a constant scattering length. Hydrogen scatters about as strongly as carbon, but with the opposite sign. Missing hydrogens are therefore added implicitly at their bond distance from the atoms that carry them.
- **The solvent.** The solvent is an H₂O/D₂O mixture, pure D₂O by default. The excluded volume and the hydration shell both take its isotope composition.
- **Hydrogen exchange.** Implicit hydrogens bound to nitrogen, oxygen and sulfur exchange with the solvent. Those bound to carbon do not.
- **Deuterated samples.** Hydrogen and deuterium atoms written explicitly in the structure file are used as given with `hydrogens --keep`. To model a deuterated protein, write its deuterium atoms into the file.

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

The solvent options must be set before the structure is loaded:

| Command line       | pyAUSAXS setting              | Description |
|--------------------|-------------------------------|-------------|
|`solv --d2o`        | `neutron_solvent_density`     | The solvent. On the command line, its D₂O fraction from 0 to 1. In Python, its scattering length density in fm/Å³: −0.0559 for H₂O and 0.6363 for D₂O, linear in between. The default is pure D₂O. |
|`solv --exchange`   | `neutron_exchange_efficiency` | The fraction of the exchangeable hydrogens that take the isotope composition of the solvent, from 0 to 1. The default is 1, full exchange. |

Apart from these, `ausaxs sans` takes the same arguments as `ausaxs fit`. The full documentation is available [here](saxs-fitter.md).
