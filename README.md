# Todo
 * [ ] Grid tests: Do all of them with all three placement strategies
 * [ ] Vector3: Template this class so it can also be used with unsigned ints. Use it instead of raw vectors in Grid methods. 
 * [ ] IntensityFitter: Locally stored fit does not contain plots. Currently fixed by returning another fit. 
 * [ ] Protein: Figure out how getMass should work in without effective charge enabled
 * [ ] EM: Figure out why the test partial_histogram_manager::comparison with standard approach doesn't work with data/A2M/A2M_na.ccp4. Probably something to do with assumed negative staining?
 * [ ] Slice: Change storage to be Matrix<T>* such that it is always kept up to date. Currently Slices are invalidated when Matrix data is changed (like with extend). 
 * [x] EMFitter: We have a lot of measurements that fit their pdb structures really well. Try to generate maps from the pdb structures and see if they also fit their maps well. 
 * [ ] EMFitter: Fit 12752, 12753 trypsin-activated
 * [ ] EMFitter: Handle absolute scale properly.
 * [x] Use nothrow
 * [ ] gentag alt for de andre filer
 * [ ] flow diagram med hvad vi rent faktisk laver
 * [ ] General: Consistency check of DrhoM
 * [ ] EM: Research and implement Electron Transfer Function (if not too difficult)
 * [ ] Atom: Add define statement controlling safety checks on set_element, get_mass and get_Z
 * [ ] Atom: Rename "effective charge" to "relative charge". Apply this change globally. 
 * [ ] EM: Do a better job of simulating experimental data (uncertainties, Gaussian noise, better spacing). Do a check on the voxel sizes and skip every Nth pixel if it is too small. 
 * [ ] IO: Support multiple terminate statements
 * [ ] Memory test all other executables.
 * [ ] Constants: Rethink how to determine charge densities for arbitrary ligands
 * [ ] Matrix & Slice: Remove debug segfaults for out-of-bounds indices

# Stuff to consider
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
 * ROOT (compile options: `cmake -DCMAKE_INSTALL_PREFIX=<install> -Dminuit2=ON -DCMAKE_CXX_STANDARD=17 -Dbuiltin_gsl=ON <source>`)
 * Boost (Very minor dependency, consider removing it entirely.)
 * Elements
 * CLI11
 * catch2 for tests

# FITTING:
FOXS SAXS fitting program
ATSAS CRYSOL

# Other personal notes
## Articles
 * Joachim (Jochen) Hub - molecular dynamics
 * Aquaporin DDM
 * B. Jacrot rep prog phys ~ 1980
 * Peter Zipper ~1980
 * ATSAS Crysol
