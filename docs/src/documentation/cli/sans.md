---
title: Fitting atomic models to SANS data - AUSAXS
description: Calculate the SANS profile of a PDB or mmCIF structure and fit it to experimental small-angle neutron scattering data with ausaxs sans, with control over the D2O fraction of the solvent and hydrogen exchange.
---

# Fitting atomic models to SANS data

`ausaxs sans` is [`ausaxs fit`](fit.md) with neutron scattering lengths instead of X-ray form factors. It takes the same input files and arguments, and writes the same output, so this page only lists the differences:  
`ausaxs sans <model file> <sans file>`

For what SANS adds over SAXS, see the [SANS application](../../applications/sans.md).

## Differences from `ausaxs fit`
- **Output folder.** The default output folder is `output/sans/` instead of `output/saxs/`.
- **Units.** Charges and densities are reported as scattering lengths in fm and fm/Å³, instead of electrons.
- **Probe.** The tool always uses neutrons, even if a settings file says otherwise.
- **Solvent.** The `solv` group has two extra arguments for the composition of the solvent, listed below.

## `solv` group
In addition to the arguments of the [`solv` group of `ausaxs fit`](fit.md#solv-group):

| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`--d2o`             | Set the D₂O fraction       | The D₂O fraction of the H₂O/D₂O solvent, from 0 to 1. This sets the scattering length density of the bulk solvent, the hydration shell, and the exchanged hydrogens. Defaults to 1, pure D₂O. |
|`--exchange`        | Set the exchange efficiency | The fraction of the exchangeable hydrogens, those bound to N, O or S, that take the isotope composition of the solvent, from 0 to 1. Defaults to 1, full exchange. |

The solvent composition used is printed at the start of every run.

As an example, `ausaxs sans model.pdb data.dat solv --d2o 0.42 --exchange 0.9` fits a measurement in 42% D₂O, where 90% of the exchangeable hydrogens have exchanged.
