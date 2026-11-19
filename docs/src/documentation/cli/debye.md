---
title: Debye scattering of atomic structures - AUSAXS
description: Calculate the Debye scattering of crystals, nanoparticles and other atomic structures without solvent with ausaxs debye, optionally as a powder pattern in 2θ or fitted to a measurement. Covers the output files and all arguments.
---

# Debye scattering of atomic structures

`ausaxs debye` calculates the scattering of a structure with the Debye equation, without any solvent: there is no excluded volume and no hydration shell. It is meant for crystals, nanoparticles and other samples measured in vacuum or air, see the [crystals and nanoparticles](../../applications/crystals.md) and [powder diffraction](../../applications/powder-diffraction.md) applications.

The tool is built on [`ausaxs fit`](fit.md), and shares many of its arguments. This page lists the arguments of its own, and refers to the `fit` page for the rest. The only required argument is the structure file:  
`ausaxs debye <structure file>`

This calculates the scattering intensity for 1000 evenly spaced $q$ values from $10^{-4}$ to 1 Å⁻¹. Unlike `ausaxs fit`, the $q$-range is not limited to 1 Å⁻¹, so wide-angle patterns can be calculated with e.g. `--qmax 10`.

If a measurement is given as well, the intensity is calculated at its $q$ values instead, and fitted to it with a scale factor and a constant background:  
`ausaxs debye <structure file> <data file>`

Any water molecules in the structure file are discarded.

## Output
By default all output is saved in the `output/debye/` folder, which can be changed with the `--output` argument.

- Without a measurement, the profile is written to **debye.dat** in `output/debye/simulated/<structure-name>/`, with the columns $q$ and $I$.
- With a measurement, the fit is written to **ausaxs.fit** in `output/debye/<data-name>/`, along with a **report.txt**. Both are described on the [`ausaxs fit` page](fit.md#output).

The header of each file states the number of atoms and how they scatter. With `--wavelength`, both files get an extra column with the scattering angle $2\theta$ in degrees.

## Arguments
Arguments are given the same way as for `ausaxs fit`, with some of them in [argument groups](fit.md#argument-groups). Run `ausaxs debye --help` for a summary.

### Base group
In addition to the general arguments of the [base group of `ausaxs fit`](fit.md#base-group), such as `--output`, `--threads` and `--gpu`:

| Argument           | Name                       | Description |
|--------------------|----------------------------|-------------|
|`--qmin`            | Set minimum q              | The lower limit of the $q$-range in Å⁻¹. Defaults to $10^{-4}$. With a measurement, it only narrows the measured range. |
|`--qmax`            | Set maximum q              | The upper limit of the $q$-range in Å⁻¹. Defaults to 1. With a measurement, it only narrows the measured range. |
|`-n`, `--count`     | Set the number of q values | The number of $q$ values, evenly spaced from `qmin` to `qmax`. Defaults to 1000. Not used with a measurement. |
|`--wavelength`      | Set the wavelength   | Add a column with the scattering angle $2\theta$ in degrees for this wavelength in Å, from $q = \frac{4\pi}{\lambda}\sin\theta$. |
|`--bin-width`       | Set the bin width          | The width of the distance bins in Å. The default should be sufficient for most cases. |
|`--points`          | Use point scatterers       | Treat every atom as a point scatterer instead of using its form factor. Each atom is weighted by its number of electrons, or its scattering length with `--neutron`, including its implicit hydrogens and occupancy. |
|`--neutron`         | Use neutrons               | Use neutron scattering lengths instead of X-ray form factors. |

### Groups
| Group              | Description |
|--------------------|-------------|
|`config`            | As for [`ausaxs fit`](fit.md#config-group). Without a measurement, a settings file is looked for in the folder of the structure file instead. |
|`molecule`          | As for [`ausaxs fit`](fit.md#molecule-group). |
|`hydrogens`         | As for [`ausaxs fit`](fit.md#hydrogens-group). |
|`data`              | Only the `--unit` and `--skip` arguments of the [`data` group of `ausaxs fit`](fit.md#data-group). The $q$-range is set with `--qmin` and `--qmax` in the base group instead. |
|`solv`              | Only the `--d2o` and `--exchange` arguments of the [`solv` group of `ausaxs sans`](sans.md#solv-group), which set the isotope composition of the exchangeable hydrogens. Only used with `--neutron`. |

As an example, `ausaxs debye Ag_crystal.cif --qmax 8.15 --wavelength 1.5406` calculates the Cu Kα powder pattern of a silver nanoparticle.
