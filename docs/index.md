---
title: "AUSAXS: Debye-equation SAXS, WAXS and SANS calculator"
description: Documentation for AUSAXS, a fast small-angle scattering calculator built on the Debye equation. Calculate SAXS, WAXS and SANS profiles from PDB, mmCIF, CIF and XYZ structures, and fit atomic models and EM maps to experimental data.
---

# AUSAXS documentation

AUSAXS implements a highly optimized, general-purpose Debye equation kernel for the CPU and the GPU. Built on this kernel, it accurately predicts small-angle X-ray scattering (SAXS), wide-angle X-ray scattering (WAXS) and small-angle neutron scattering (SANS) profiles of proteins and other biomolecules, and fits them to experimental small-angle scattering data. The same kernel computes the Debye scattering of crystals, nanoparticles and arbitrary point scatterers, with or without atomic form factors. Structures are read from PDB, mmCIF, CIF and XYZ files. AUSAXS runs on Linux, macOS and Windows, and is available both as a command-line program and as a Python package:
```bash
pip install pyausaxs
```

## User guides
All tools are run through the single `ausaxs` executable as `ausaxs <tool> [options]`.

- **[SAXS fitting](saxs-fitter.md)** (`ausaxs fit`): Calculate the SAXS profile of a protein or other macromolecule from its atomic model, and test whether it is consistent with experimental SAXS data. The hydration shell and excluded volume are modelled explicitly from the atomic coordinates.
- **[EM map validation](em-validation.md)** (`ausaxs em`): Validate cryo-EM maps against solution SAXS data using dummy structures built from the map itself.
- **[Rigid-body refinement](rigidbody.md)** (`ausaxs rigidbody`): Refine multi-domain proteins and complexes against SAXS data, driven by configuration scripts. This tool is still in active development.
- **[Plotting scripts](plotting-scripts.md)**: Turn the output of the command-line tools into figures.

All three tools are also available through a graphical interface and a Python API in the [pyAUSAXS](https://github.com/AUSAXS/pyAUSAXS) package. Its [examples](https://github.com/AUSAXS/pyAUSAXS/tree/master/examples) cover everything from a basic Debye calculation to model fitting and crystal scattering.

## Getting started
See the [installation guide](installation.md) to install AUSAXS with `pip`, download a precompiled binary, or compile it from source. If you have issues running the programs, please take a look at the [troubleshooting](troubleshooting.md) page.

## Citing AUSAXS
If you use AUSAXS in your work, please cite:

- K. Lytje and J. S. Pedersen. Small-angle X-ray scattering profile calculation for high-resolution models of biomacromolecules. _J. Appl. Cryst._ (2025). doi: [10.1107/S160057672500562X](https://doi.org/10.1107/S160057672500562X)

For EM map validation, please also cite:

- Validation of electron-microscopy maps using solution small-angle X-ray scattering. _Acta Cryst. D_ (2024). doi: [10.1107/S2059798324005497](https://doi.org/10.1107/S2059798324005497)
