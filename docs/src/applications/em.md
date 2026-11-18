---
title: Validating cryo-EM maps with SAXS - AUSAXS
description: How AUSAXS validates cryo-EM maps against solution SAXS data, by fitting hydrated dummy models of the map at every threshold level and mapping out the chi-square landscape.
journey: apps
hide:
  - toc
---

# EM map validation

<span class="ausaxs-name"><span class="au">AU</span>SAXS</span> validates a cryo-EM map against solution SAXS data. It tells you whether the map agrees with the structure of the molecule in solution, and at which threshold level.

## Why combine EM and SAXS
- **The threshold.** A cryo-EM map has no absolute density scale. The shape of the molecule depends on the threshold level at which the map is cut, and the choice of level is largely a matter of judgement.
- **The solution state.** A cryo-EM map is reconstructed from molecules frozen on a grid, and averaged over particles that were selected and classified. Solution SAXS measures every molecule in the sample, in its native buffer.
- **The comparison.** Cutting the map at a threshold level gives a shape whose SAXS profile can be calculated and compared with the data. Scanning over the levels gives the χ² as a function of the level. A clear minimum means that the map agrees with the molecule in solution at that level; no good fit at any level points to a difference between the two states.

The method is described in [Lytje & Pedersen (2024)](https://doi.org/10.1107/S2059798324005497).

## Example
Here we validate the cryo-EM map of the SARS-CoV-2 main protease, [EMD-24889](https://www.ebi.ac.uk/emdb/EMD-24889), [emd_24889.map](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/emd_24889.map), against a SAXS measurement of it, [SASDJG5.dat](https://raw.githubusercontent.com/AUSAXS/AUSAXS/master/tests/files/SASDJG5.dat), and plot the result:

```bash
ausaxs em emd_24889.map SASDJG5.dat
ausaxs plot --folder output/em_fitter/SASDJG5/emd_24889
```

![The chi-square of the map falls steeply with the threshold level, reaching a broad minimum near 7.4 sigma; the fit at the minimum follows the measured SAXS data.](../assets/figures/em-validation.svg#only-light){ width="100%" }
![The chi-square of the map falls steeply with the threshold level, reaching a broad minimum near 7.4 sigma; the fit at the minimum follows the measured SAXS data.](../assets/figures/em-validation-dark.svg#only-dark){ width="100%" }

For a worked example on the same data, see the [EM validation tutorial](https://sastutorials.org/SAXS_EM/SAXS_EM.html). The full documentation for the `em` program is available [here](../documentation/cli/em-validation.md).
