# Fitting atomic models to SANS data

`ausaxs sans` is [`ausaxs fit`](fit.md) with neutron scattering lengths instead of X-ray form factors. It takes the same input files and arguments, and writes the same output, so this page only lists the differences:  
`ausaxs sans <model file> <sans file>`

For what SANS adds over SAXS, see the [SANS application](https://ausaxs.github.io/applications/sans/).

## Differences from `ausaxs fit`
- **Output folder.** The default output folder is `output/sans/` instead of `output/saxs/`.
- **Units.** Charges and densities are reported as scattering lengths in fm and fm/Å³, instead of electrons.
- **Probe.** The tool always uses neutrons, even if a settings file says otherwise.
- **Solvent.** The `solv` group has two extra arguments for the composition of the solvent, listed below.
- **Hydrogens.** `hydrogens --keep` is the default: the hydrogens in the file are used as they are, since neutrons resolve their positions. A file without any gets implicit hydrogens instead, with a warning, as they make the scattering a few percent too low at q ≳ 0.2 Å⁻¹. Use `hydrogens --discard` to always use implicit hydrogens.

## `solv` group
In addition to the arguments of the [`solv` group of `ausaxs fit`](fit.md#solv-group):

| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`--d2o`             | Set the D₂O fraction       | The D₂O fraction of the H₂O/D₂O solvent, from 0 to 1. This sets the scattering length density of the bulk solvent, the hydration shell, and the exchanged hydrogens. Defaults to 1, pure D₂O. |
|`--exchange`        | Set the exchange efficiency | The fraction of the exchangeable implicit hydrogens, those bound to N, O or S, that take the isotope composition of the solvent, from 0 to 1. |

> [!WARNING]
> `--d2o` defaults to pure D₂O. For a measurement in H₂O or in a mixture, always give its D₂O fraction.

The solvent composition used is printed at the start of every run.

`--d2o` and `--exchange` only change the implicit hydrogens that AUSAXS adds itself. Explicit hydrogens are used exactly as given in the file, so for a measurement in D₂O, the exchangeable sites must already be deuterium there.

As an example, `ausaxs sans model.pdb data.dat solv --d2o 0.42 --exchange 0.9` fits a measurement in 42% D₂O, where 90% of the exchangeable hydrogens have exchanged.
