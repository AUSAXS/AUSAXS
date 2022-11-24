# Todo
 * [ ] CCP4 reader: Major bug discovered. Ensure reading is *always* done along same axes as specified in the header. 
 * [ ] DebyeLookupTable: Non-default tables seems to be broken. Disable default tables to debug. 
 * [ ] RhoM: Calculate total weight of protein, multiply (divide) by Avogadro constant, divide by volume to get average density. Should be 1.4 or something like that. 
 * [ ] EM: R factors http://pd.chem.ucl.ac.uk/pdnn/refine1/rfacs.htm
 * [ ] EM: Compare maps from different simulation methods
 * [ ] ResidueParser: Consider creating 1 parsed file per residue and then dynamically loading only those that are needed for the current file. 
 * [ ] Protein: Figure out how getMass should work in without effective charge enabled
 * [ ] EM: Figure out why the test partial_histogram_manager::comparison with standard approach doesn't work with data/A2M/A2M_ma.ccp4. Probably something to do with assumed negative staining?
 * [ ] Slice: Change storage to be Matrix<T>* such that it is always kept up to date. Currently Slices are invalidated when Matrix data is changed (like with extend). 
 * [ ] EMFitter: Handle absolute scale properly.
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
 * The grid map stores chars = 1 byte or 8 bits, but we only use 5 different values (0, 'h', 'H', 'a', 'A'). If we can reduce this to 4 values, it can be stores in just 2 bits, in which case the remaining 6 bits can be used by other threads. That's literally free multithreading. Possible race condition. 

## ScatteringHistogram:
 * Optional argument of q-values to calculate I(q) for - this would remove the necessity of splicing in the IntensityFitter
 * Take a closer look at the form factor
 * Convert `a` to a more sensible output (units)

## Protein: 
 * When calculating the volume of the acids, the calculation is simply delegated to each individual body. However, if an amino acid happens to be cut in two halves in two different bodies, its volume is counted twice. 

# Dependencies
Maybe bundle them somehow to make it easier to install?
 * Elements
 * CURL
 * OpenSSL (CURL)
 * binutils-dev (backward-cpp stacktraces)

# Other personal notes
## Articles
 * Joachim (Jochen) Hub - molecular dynamics
 * Aquaporin DDM
 * B. Jacrot rep prog phys ~ 1980
 * Peter Zipper ~1980
 * ATSAS Crysol
