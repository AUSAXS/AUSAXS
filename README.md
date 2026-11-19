<h1>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/AUSAXS/AUSAXS/media/title_light.png">
  <img alt="AUSAXS: small-angle scattering calculator" src="https://raw.githubusercontent.com/AUSAXS/AUSAXS/media/title_dark.png">
</picture>
</h1>

[![release](https://img.shields.io/github/v/release/AUSAXS/AUSAXS)](https://github.com/AUSAXS/AUSAXS/releases/latest)
[![PyPI](https://img.shields.io/pypi/v/pyausaxs)](https://pypi.org/project/pyausaxs/)
[![coverity](https://scan.coverity.com/projects/30350/badge.svg)](https://scan.coverity.com/projects/ausaxs-ausaxs)
[![license](https://img.shields.io/github/license/AUSAXS/AUSAXS)](LICENSE)

**AUSAXS** implements a highly optimized, general-purpose Debye equation kernel for both the CPU and GPU. Built on this kernel, it accurately predicts small-angle X-ray scattering (SAXS), wide-angle X-ray scattering (WAXS) and small-angle neutron scattering (SANS) profiles of proteins and other biomolecules, and fits them to experimental small-angle scattering data. The same kernel computes the Debye scattering of crystals, nanoparticles and arbitrary point scatterers, with or without atomic form factors. AUSAXS runs on Linux, macOS and Windows, and is available both as a command-line program and as a Python package:
```bash
python -m pip install pyausaxs
```

# Main features
- **[SAXS profile calculation and fitting](https://ausaxs.github.io/AUSAXS/documentation/cli/fit/)**: Calculate the SAXS curve of a structure and test whether it is consistent with experimental data, adjusting only solvent and scaling parameters. The hydration shell and excluded volume are modelled explicitly from the atomic coordinates.
- **[SANS profile calculation and fitting](https://ausaxs.github.io/AUSAXS/documentation/cli/sans/)**: The same calculation for small-angle neutron scattering, with neutron scattering lengths, hydrogen/deuterium exchange and D₂O contrast variation.
- **[EM map validation](https://ausaxs.github.io/AUSAXS/documentation/cli/em-validation/)**: Compare cryo-EM maps (MRC/CCP4) against SAXS data using dummy structures built from the map itself.
- **[Rigid-body refinement](https://ausaxs.github.io/AUSAXS/documentation/cli/rigidbody/)** _(in development)_: Refine multi-domain proteins and complexes against SAXS data, exploiting symmetry relations to accelerate the calculation. Driven by a composable configuration language for fine-grained control over parameters, loops, and per-step output.
- **Debye scattering of crystals and point scatterers**: Evaluate the Debye equation at any _q_ for crystals, nanoparticles and other materials read from CIF or XYZ files, or for arbitrary point scatterers given as coordinates and weights, with or without atomic form factors. The gradient is also available, for use in optimization and machine-learning workflows. Available through the [Python package](https://github.com/AUSAXS/pyAUSAXS).
- **GPU acceleration**: Optional GPU backend for NVIDIA and AMD cards, installed and checked with `ausaxs gpu`.

# Performance
AUSAXS is the fastest SAXS program by a wide margin. Even systems containing millions of atoms can be evaluated in seconds on consumer GPUs.

The GPU is optional, though. On the CPU alone, AUSAXS still outpaces every other program: a good consumer CPU handles about 300 000 atoms in 10 seconds.

The benchmarks, against both SAXS programs and Debye calculators, are in the [documentation](https://ausaxs.github.io/AUSAXS/why-ausaxs/#speed).

# Documentation
Each program has a detailed user guide in the [documentation](https://ausaxs.github.io/AUSAXS/). Examples of the Python interface, from a basic Debye calculation to model fitting and crystal scattering, are in the [pyAUSAXS examples](https://github.com/AUSAXS/pyAUSAXS/tree/master/examples).

# Installation
**Python (recommended)**: The easiest way to get started is the Python package, which installs the `ausaxs` command-line tool in one step:
```bash
python -m pip install pyausaxs
```
No need to track down a downloaded binary: just `pip install` and you're ready to go. The package also exposes most features programmatically, making AUSAXS easy to integrate into external workflows. See the [pyAUSAXS](https://github.com/AUSAXS/pyAUSAXS) repository for details.

**Binaries**: Precompiled binaries are also available directly in the [releases](https://github.com/AUSAXS/AUSAXS/releases).

**From source**: Follow the [compilation guide](https://ausaxs.github.io/AUSAXS/getting-started/compiling/) to build it yourself in a few steps.

# Contributing
Encountering problems, have feedback or suggestions, or considering contributing? Please check out the [contributor guidelines](CONTRIBUTING.md).

# Citing AUSAXS
If you use AUSAXS in your work, please cite:
* K. Lytje and J. S. Pedersen. Small-angle X-ray scattering profile calculation for high-resolution models of biomacromolecules. _J. Appl. Cryst._ (2025)  
(doi: [10.1107/S160057672500562X](https://doi.org/10.1107/S160057672500562X))

For EM map validation, please also cite:
* Validation of electron-microscopy maps using solution small-angle X-ray scattering. _Acta Cryst. D_ (2024)  
(doi: [10.1107/S2059798324005497](https://doi.org/10.1107/S2059798324005497))

_This project is licenced under the GNU Lesser General Public Licence v3.0. Supported by grant 1026-00209B from the Independent Research Fund Denmark._
