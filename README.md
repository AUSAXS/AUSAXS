![title_light](../media/title_dark.png?raw=true#gh-light-mode-only)
![title_dark](../media/title_light.png?raw=true#gh-dark-mode-only)

# Main features
- **Simple foundation**: We have implemented the methods in the simplest possible way, making as few assumptions about your data as possible. With the Debye equation as the basis for the scattering profiles, the only loss of accuracy is through the histogram approximation, where we support using both weighted and unweighted bins depending on your preferences. By implementing the technique in modern C++ with efficiency in mind, we have managed to achieve some of the [best performance available](https://github.com/klytje/AUSAXS/blob/media/benchmark.png).
- **Fitting of high-resolution models to SAXS curves**: Fit atomic structure files using experimental SAXS data using an efficient implementation of the Debye equation. Various options are available regarding the handling of both the excluded volume and hydration shell. 
- **Validation of electron microscopy maps**: Validate EM maps using experimental SAXS data. By using the information contained within the EM map itself, dummy structures can be constructed and compared against the SAXS data. Though various other implementations doing something similar are already available ([scipion](scipion.i2pc.es), [denss](https://tdgrant.com/)), ours is the only one that manages to consistenly achieve single-digit $\chi^2$ values for matching experimental datasets. 
- **Rigidbody optimization**: _(Still under development)_ Perform self-consistent and customizable rigidbody optimizations, generating a new hydration shell for each step. Optional calibration with scattering curves predicted by molecular dynamics simulations can limit the number of free parameters to just 2, dramatically reducing the capability of overfitting.

User-guides to all of these programs can be found in the [wiki](https://github.com/klytje/AUSAXS/wiki).

# Installation
## Download precompiled binaries
The fastest way to get started is using the most recent precompiled executables available in the [releases](https://github.com/klytje/AUSAXS/releases). Alternatively you can follow the next section to compile the library yourself. 

## Compile from source
The software can easily be compiled from source with only a few steps. GCC v11+, Clang v15+, and MSVC 2022+ are supported, though GCC is the preferred option for optimal efficiency.

### Linux
1. Make sure you have the prerequisites installed  
`apt-get install cmake make g++ libcurl4-openssl-dev`

2. Clone this repository  
`git clone https://github.com/klytje/AUSAXS.git`.

3. Run the build command  
`make build`

4. Compile your choice of executable  
`make intensity_fitter`

### Windows
1. Make sure CURL and OpenSSL are available on your system, e.g. through vcpkg

2. Download or clone this repository
`git clone https://github.com/klytje/AUSAXS.git`.

3. Compile your choice of executable. Note that this is very memory-intensive with the MSVC compiler, requiring 12GB+ of available memory due to their inefficient handling of constant expressions. 

### Mac
Note that Mac is not officially supported since I don't have such a machine available for testing. Make sure you have the most recent Clang compiler and `curl` available on your system, and then follow the Linux steps. 

# Dependencies
Manual dependencies:
*	[CURL](https://github.com/curl/curl), for dynamically downloading configuration files for unknown complexes and ligands. 
	*	[OpenSSL](https://github.com/openssl/openssl)

The following are automatically fetched by CMake:
*	[Elements](https://github.com/cycfi/elements), for the graphical user interface.
	*	[Native File Dialog Extended](https://github.com/btzy/nativefiledialog-extended), for the file dialogues. 
*	[Generalized Constant Expression Math](https://github.com/kthohr/gcem), for constexpr math support beyond GCC. 
*	[BSThreadPool](https://github.com/bshoshany/thread-pool), for multithreading. 
*	[CLI11](https://github.com/CLIUtils/CLI11), for interpreting command line arguments. 
*	[dlib](https://github.com/davisking/dlib), for their fitting routine. 
*	[Catch2](https://github.com/catchorg/Catch2), for running tests.
*	[backward-cpp](https://github.com/bombela/backward-cpp), for better stacktraces in debug builds. 

# Disclaimer
This project was supported by grant 1026-00209B from the Independent Research Fund Denmark. 
