![title_light](../media/title_dark.png?raw=true#gh-light-mode-only)
![title_dark](../media/title_light.png?raw=true#gh-dark-mode-only)

[![release](https://img.shields.io/github/v/release/AUSAXS/AUSAXS)](https://github.com/AUSAXS/AUSAXS/releases/latest)
[![coverity](https://scan.coverity.com/projects/30350/badge.svg)](https://scan.coverity.com/projects/ausaxs-ausaxs)
[![license](https://img.shields.io/github/license/AUSAXS/AUSAXS)](LICENSE)

**AUSAXS** is a modern C++20 small-angle X-ray scattering (SAXS) calculator built on the Debye equation, written to be modular and easily extendable. It is highly efficient (see the [benchmark](https://github.com/AUSAXS/AUSAXS/blob/media/benchmark.png)) and available on Linux, macOS, and Windows.

Thanks to its _weighted-bin_ implementation of the Debye equation, AUSAXS can accurately predict the pure Debye scattering in the wide-angle (WAXS) regime, up to arbitrary _q_. At smaller angles, our novel explicit-solvent hydration-shell and grid-based excluded-volume modeling enable reliable validation of structures against SAXS data. 

# Main features
- **[Structural validation](https://github.com/AUSAXS/AUSAXS/wiki/SAXS-fitter)**: Test whether an atomic model is consistent with experimental SAXS data, adjusting only solvent and scaling parameters.
- **[EM map validation](https://github.com/AUSAXS/AUSAXS/wiki/Validation-of-EM-maps)**: Compare electron-microscopy maps against SAXS data using dummy structures built from the map itself.
- **[Rigidbody optimization](https://github.com/AUSAXS/AUSAXS/wiki/Rigidbody-optimization)** _(in development)_: Refine multi-body assemblies against SAXS data, exploiting symmetry relations to accelerate the calculation. Driven by a composable configuration language for fine-grained control over parameters, loops, and per-step output.

Each program has a detailed user guide in the [wiki](https://github.com/AUSAXS/AUSAXS/wiki).

# Installation
**Python (recommended)**: The easiest way to get started is the Python wrapper, which installs the `ausaxs` command-line tool in one step:
```bash
pip install pyausaxs
```
No need to track down a downloaded binary — just `pip install` and you're ready to go. The wrapper also exposes most main features programmatically, making AUSAXS easy to integrate into external workflows. See the [pyausaxs](https://github.com/AUSAXS/pyAUSAXS) repository for details.

**Binaries**: Precompiled binaries are also available directly in the [releases](https://github.com/AUSAXS/AUSAXS/releases).

**From source**: Follow the [compilation guide](https://github.com/AUSAXS/AUSAXS/wiki/Compilation-&-installation) to build it yourself in a few steps.

# Contributing
Encountering problems, have feedback or suggestions, or considering contributing? Please check out the [contributor guidelines](CONTRIBUTING.md).

# References
* Small-angle X-ray scattering profile calculation for high-resolution models of biomacromolecules  
(doi: [10.1107/S160057672500562X](https://doi.org/10.1107/S160057672500562X))
* Validation of electron-microscopy maps using solution small-angle X-ray scattering  
(doi: [10.1107/S2059798324005497](https://doi.org/10.1107/S2059798324005497))

_This project is licenced under the GNU Lesser General Public Licence v3.0. Supported by grant 1026-00209B from the Independent Research Fund Denmark._
