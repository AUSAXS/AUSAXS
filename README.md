![title_light](../media/title_dark.png?raw=true#gh-light-mode-only)
![title_dark](../media/title_light.png?raw=true#gh-dark-mode-only)

[![release](https://img.shields.io/github/v/release/AUSAXS/AUSAXS)](https://github.com/AUSAXS/AUSAXS/releases/latest)
[![coverity](https://scan.coverity.com/projects/30350/badge.svg)](https://scan.coverity.com/projects/ausaxs-ausaxs)

# Main features 
- **Simple foundation**: We have implemented the methods in the simplest possible way, making as few assumptions about your data as possible. With the Debye equation as the basis for the scattering profiles, the only loss of accuracy is through the histogram approximation, where we support using both weighted and unweighted bins depending on your preferences. By implementing the technique in modern C++ with efficiency in mind, we have managed to achieve some of the [best performance available](https://github.com/klytje/AUSAXS/blob/media/benchmark.png).
- **Fitting of high-resolution models to SAXS curves**: Fit atomic structure files using experimental SAXS data using an efficient implementation of the Debye equation. Various options are available regarding the handling of both the excluded volume and hydration shell. 
- **Validation of electron microscopy maps**: Validate EM maps using experimental SAXS data. By using the information contained within the EM map itself, dummy structures can be constructed and compared against the SAXS data, serving as a quick quality check on the conformation of the map. 
- **Rigidbody optimization**: _(Still under development)_ Perform self-consistent and customizable rigidbody optimizations, generating a new hydration shell for each step. Optional calibration with scattering curves predicted by molecular dynamics simulations can limit the number of free parameters to just 2, dramatically reducing the capability of overfitting.

User-guides to all of these programs can be found in the [wiki](https://github.com/klytje/AUSAXS/wiki).

# Installation
The fastest way to get started is using the most recent precompiled binaries available in the [releases](https://github.com/klytje/AUSAXS/releases).  

Alternatively you can follow the [compilation guide](https://github.com/AUSAXS/AUSAXS/wiki/Compilation-&-installation) to compile it yourself in just a few simple steps. 

# Contributing
Are you encountering problems, have feedback or suggestions, or are you considering contributing to the project? Please check out the [contributor guidelines](CONTRIBUTING.md).  

# References
* Small-angle X-ray scattering profile calculation for high-resolution models of biomacromolecules  
(doi: [10.1107/S160057672500562X](https://doi.org/10.1107/S160057672500562X))
* Validation of electron-microscopy maps using solution small-angle X-ray scattering (doi: [10.1107/S2059798324005497](https://doi.org/10.1107/S2059798324005497))

Another article documenting the rigidbody optimizer is currently under development. 

_This project is licenced under the GNU Lesser General Public Licence v3.0. Alternative licencing arrangements can be discussed upon request. Supported by grant 1026-00209B from the Independent Research Fund Denmark._
