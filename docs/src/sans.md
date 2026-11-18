---
title: SANS profiles of biomolecules - AUSAXS
description: How AUSAXS calculates small-angle neutron scattering (SANS) profiles of proteins from atomic models, with neutron scattering lengths, H2O/D2O solvent contrast and hydrogen-deuterium exchange.
journey: apps
---

# Small-angle neutron scattering (SANS)

AUSAXS calculates SANS profiles with the same pipeline as [SAXS](saxs.md): the atoms, the excluded volume and the hydration shell are modelled in the same way, and the X-ray form factors are replaced by neutron scattering lengths.

## What changes for neutrons
- **Scattering lengths.** A nucleus scatters as a point, so a single atom has a constant scattering length. Hydrogen has a negative scattering length, and the hydrogens bound to each heavy atom are placed at their real bond distance. This gives groups such as CH₂ a strong $q$-dependence.
- **Solvent contrast.** The solvent is an H₂O/D₂O mixture, set by its scattering length density. The default is pure D₂O. The hydration shell takes the isotope composition of the solvent.
- **Hydrogen–deuterium exchange.** Implicit hydrogens bound to nitrogen, oxygen and sulfur exchange with the solvent. The exchange efficiency is a setting, with full exchange as the default. Hydrogens bound to carbon do not exchange.
- **Deuterated samples.** Hydrogen and deuterium atoms written explicitly in the structure file are used as given. To model a deuterated protein, write its deuterium atoms into the file.

## Running it
From the command line, `ausaxs sans` takes the same arguments as [`ausaxs fit`](saxs-fitter.md), plus two for the solvent:
```bash
ausaxs sans model.pdb data.dat solv --d2o 0.42
```

| Argument           | Description |
|--------------------|-------------|
|`solv --d2o`        | The D₂O fraction of the solvent, from 0 to 1. The default is 1, pure D₂O. |
|`solv --exchange`   | The fraction of the exchangeable hydrogens that take the isotope composition of the solvent, from 0 to 1. The default is 1, full exchange. |

From Python, select neutrons before loading the structure:
```python
import pyausaxs as ausaxs
ausaxs.settings.set("radiation", "neutron")
mol = ausaxs.create_molecule("model.pdb")
mol.hydrate()
result = mol.fit("data.dat")
```
