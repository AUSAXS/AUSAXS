![title_light](../media/title_dark.png?raw=true#gh-light-mode-only)
![title_dark](../media/title_light.png?raw=true#gh-dark-mode-only)

step 1: udskift formfaktor i oprindeligt program
step 2: kig på bidraget af hydreringen til intensiteten

# Main features
- **Simple foundation**: We strive to use only simple methods and techniques, making as few assumptions as possible. By implementing the methods in modern C++ with efficiency in mind, we have managed to achieve some of the best performance available.
- **Rigidbody optimization**: Perform self-consistent and customizable rigidbody optimizations, generating a new hydration shell for each step. Optional calibration with scattering curves predicted by molecular dynamics simulations can limit the number of free parameters to just 2, dramatically reducing the capability of overfitting.   
- **Validation of electron microscopy maps**: Validate EM maps using experimental SAXS data. 
- **Fitting of high-resolution models to SAXS curves**: Fit atomic structure files using experimental SAXS data. 

# Installation
## Compile from source
The software can easily be compiled from source with only a few steps.

### Linux
1. Make sure you have the prerequisites installed  
`apt-get install cmake make g++ libcurl4-openssl-dev`

2. Clone this repository  
`git clone https://github.com/klytje/SAXS.git`.

3. Run the build command  
`make build`

4. Compile your choice of executable  
`make intensity_fitter`

## Download precompiled binaries
Precompiled binaries are available **here**. 

# Design philosophy

# Todo
* 	[ ] HistogramManager: Consider calculating d^2 and then taking the square root afterwards. This should be a significant speedup. 
* 	[ ] Consider implementing a unit class to avoid unit conversion errors in functions
* 	[ ] Redo the storage of Image to be a single 3D vector owned by the ImageStack. 
* 	[ ] Consider implementing form factors for the most common atoms. This will also allow fitting the excluded volume to the MD simulations (probably), thus avoiding the background fitting parameter. 
*	[ ] Matrix: Change vector of vectors constructor to be row-majority, and use the variadic constructor for columns instead
*	[ ] EM: R factors http://pd.chem.ucl.ac.uk/pdnn/refine1/rfacs.htm
*	[ ] Protein: Figure out how getMass should work in without effective charge enabled
*	[ ] Slice: Change storage to be Matrix<T>* such that it is always kept up to date. Currently Slices are invalidated when Matrix data is changed (like with extend). 
*	[ ] General: Consistency check of DrhoM
*	[ ] IO: Support multiple terminate statements
*	[ ] Memory test all other executables.

# Stuff to consider
## RigidBody
* 	Use a combined translation + rotation matrix when rotating about a constraint. Requires the introduction of 4D vectors.

## EM
*	Consider the effects of discretization. 
*	e2pdb2mrc.py: 
 	eval "$(/home/au561871/tools/eman2/bin/conda shell.bash hook)"
	python3 ~/tools/eman2/bin/e2pdb2mrc.py
	
## Compiler flags:
*	fno-finite-math-only: Can probably be removed, not sure of performance benefits. Note that its removal would specifically break limit handling of plots where std::isinf checks are used. 

## Grid:
*	Consider how to improve culling method
*	Consider removing all bounds checks (or maybe use a compile-flag to enable them)

## ScatteringHistogram:
*	Take a closer look at the form factor
*	Convert `a` to a more sensible output (units)

## Protein: 
*	When calculating the volume of the acids, the calculation is simply delegated to each individual body. However, if an amino acid happens to be cut in two halves in two different bodies, its volume is counted twice. 

# Dependencies
*	Elements (currently manual install)
*	CURL
*	OpenSSL (CURL)
*	backward-cpp (can be removed in production)
	*	binutils-dev
*	BSThreadPool (currently manual install)
*	CLI11
*	dlib
*	Catch2

# Other personal notes
## Articles
*	Joachim (Jochen) Hub - molecular dynamics
*	Aquaporin DDM
*	B. Jacrot rep prog phys ~ 1980
*	Peter Zipper ~1980
*	ATSAS Crysol
