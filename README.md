# Todo
 * [ ] BSThreadPool: Streamline download and build process. Currently it is all manual. 
 * [ ] DebyeLookupTable: Non-default tables seems to be broken. Disable default tables to debug. 
 * [x] RhoM: Calculate total weight of protein, multiply (divide) by Avogadro constant, divide by volume to get average density. Should be 1.4 or something like that. 
 * [ ] EM: R factors http://pd.chem.ucl.ac.uk/pdnn/refine1/rfacs.htm
 * [ ] Protein: Figure out how getMass should work in without effective charge enabled
 * [ ] EM: Figure out why the test partial_histogram_manager::comparison with standard approach doesn't work with data/A2M/A2M_ma.ccp4. Probably something to do with assumed negative staining?
 * [ ] Slice: Change storage to be Matrix<T>* such that it is always kept up to date. Currently Slices are invalidated when Matrix data is changed (like with extend). 
 * [x] EMFitter: Handle absolute scale properly. (Doesn't matter. It is scaled to fit SAXS data anyways)
 * [ ] gentag alt for de andre filer
 * [ ] General: Consistency check of DrhoM
 * [ ] IO: Support multiple terminate statements
 * [ ] Memory test all other executables.

# Stuff to consider
## EM
 * Consider the effects of discretization. 
 * e2pdb2mrc.py: 
 	eval "$(/home/au561871/tools/eman2/bin/conda shell.bash hook)"
	python3 ~/tools/eman2/bin/e2pdb2mrc.py
	
## Compiler flags:
 * fno-finite-math-only: Can probably be removed, not sure of performance benefits. Note that its removal would specifically break limit handling of plots where std::isinf checks are used. 

## Grid:
 * Consider how to improve culling method
 * Consider removing all bounds checks (or maybe use a compile-flag to enable them)

## ScatteringHistogram:
 * Optional argument of q-values to calculate I(q) for - this would remove the necessity of splicing in the IntensityFitter
 * Take a closer look at the form factor
 * Convert `a` to a more sensible output (units)

## Protein: 
 * When calculating the volume of the acids, the calculation is simply delegated to each individual body. However, if an amino acid happens to be cut in two halves in two different bodies, its volume is counted twice. 

# Dependencies
 * Elements (currently manual install)
 * CURL
 * OpenSSL (CURL)
 * backward-cpp (can be removed in production)
	* binutils-dev
 * BSThreadPool (currently manual install)
 * CLI11
 * dlib
 * Catch2

# Other personal notes
## Articles
 * Joachim (Jochen) Hub - molecular dynamics
 * Aquaporin DDM
 * B. Jacrot rep prog phys ~ 1980
 * Peter Zipper ~1980
 * ATSAS Crysol
