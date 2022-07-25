# Todo
 * [ ] EM: Figure out why the test partial_histogram_manager::comparison with standard approach doesn't work with data/A2M/A2M_na.ccp4. Probably something to do with assumed negative staining?
 * [ ] Slice: Change storage to be Matrix<T>* such that it is always kept up to date. Currently Slices are invalidated when Matrix data is changed (like with extend). 
 * [ ] Dataset & Histogram: Basically same class, can maybe be combined somehow?
 * [ ] EMFitter: We have a lot of measurements that fit their pdb structures really well. Try to generate maps from the pdb structures and see if they also fit their maps well. 
 * [ ] EMFitter: Fit 12752, 12753 trypsin-activated
 * [ ] EMFitter: Handle absolute scale properly.
 * [ ] Use nothrow
 * [ ] gentag alt for de andre filer
 * [ ] flow diagram med hvad vi rent faktisk laver
 * [ ] General: Consistency check of DrhoM
 * [ ] EM: Research and implement Electron Transfer Function (if not too difficult)
 * [ ] Atom: Add define statement controlling safety checks on set_element, get_mass and get_Z
 * [ ] Atom: Rename "effective charge" to "relative charge". Apply this change globally. 
 * [ ] EM: Do a better job of simulating experimental data (uncertainties, Gaussian noise, better spacing). Do a check on the voxel sizes and skip every Nth pixel if it is too small. 
 * [ ] IO: Support multiple terminate statements
 * [ ] Memory test all other executables.
 * [ ] General: Determine where hydration_atoms should be stored. The Protein class seems like the best choice. 
 * [ ] Atom: Const uid 
 * [ ] Constants: Rethink how to determine charge densities for arbitrary ligands

# Stuff to consider
## Compiler flags:
 * fno-finite-math-only: Can probably be removed, not sure of performance benefits. Note that its removal would specifically break limit handling of plots where std::isinf checks are used. 

## Grid:
 * Consider simply calculating and using the bounding box at initialization as the entire grid instead of wasting memory
 * Consider how to improve culling method
 * Consider removing all bounds checks

## ScatteringHistogram:
 * Optional argument of q-values to calculate I(q) for - this would remove the necessity of splicing in the IntensityFitter
 * Take a closer look at the form factor
 * Convert `a` to a more sensible output (units)

## Body:
 * Consider removing hydration_atoms, they're not supposed to be used anyway

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
