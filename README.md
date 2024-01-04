![title_light](../media/title_dark.png?raw=true#gh-light-mode-only)
![title_dark](../media/title_light.png?raw=true#gh-dark-mode-only)

# Main features
- **Simple foundation**: We strive to use only simple methods and techniques, making as few assumptions as possible. By implementing the methods in modern C++ with efficiency in mind, we have managed to achieve some of the best performance available.
- **Rigidbody optimization**: Perform self-consistent and customizable rigidbody optimizations, generating a new hydration shell for each step. Optional calibration with scattering curves predicted by molecular dynamics simulations can limit the number of free parameters to just 2, dramatically reducing the capability of overfitting.   
- **Validation of electron microscopy maps**: Validate EM maps using experimental SAXS data. 
- **Fitting of high-resolution models to SAXS curves**: Fit atomic structure files using experimental SAXS data. 

# Installation
## Compile from source
The software can easily be compiled from source with only a few steps. GCC, Clang v15+, and MSVC are supported, though GCC is the preferred option for optimal efficiency.

### Linux
1. Make sure you have the prerequisites installed  
`apt-get install cmake make g++ libcurl4-openssl-dev`

2. Clone this repository  
`git clone https://github.com/klytje/AUSAXS.git`.

3. Run the build command  
`make build`

4. Compile your choice of executable  
`make intensity_fitter`

## Download precompiled binaries
Precompiled binaries are available **here**. 

### Windows
1. Make sure CURL and OpenSSL are available on your system, e.g. through vcpkg

2. Download or clone this repository
`git clone https://github.com/klytje/AUSAXS.git`.

3. Compile your choice of executable. Note that this is very memory-intensive with the MSVC compiler, requiring 12GB+ of available memory. 

# Design philosophy

# Dependencies
Manual dependencies:
*	CURL
*	OpenSSL

The following are automatically fetched by CMake:
*	Elements
	*	Native File Dialog Extended
*	Generalized Constant Expression Math
*	BSThreadPool
*	CLI11
*	dlib
*	Catch2 (testing only)
*	backward-cpp (debug only)
	*	binutils-dev
