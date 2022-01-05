# Todo
 * [ ] Atom: Const uid 
 * [ ] Constants: Rethink how to determine charge densities for arbitrary ligands
 * [x] Grid: Fix AxesPlacement - exact water locations, iterate over grid
 * [ ] General: Use std::array wherever possible instead of vectors. Less overhead. 
 * [ ] PartialHistogramManager: Consider redoing the data structure of CompactCoordinates (arrays?)
 * [ ] Try to derive an analytical solution of the chi2 problem. Differentiate chi2 with respect to each variable, and set each expression equal to zero. As long as it's not an iterative equation, it should be good. 
 * [ ] Make a superclass for ScatteringHistogram (maybe just use PartialHistogram?) which contains only p_tot, and defines all operations which uses only this. Then change the return type of PartialHistogramManager to this - right now some operations are ill-defined on it. 
 * [ ] Protein: something looks weird with the effective charges - run test_dist 1 & 2 from protein.cpp and look at the output
 * [ ] Grid: Something weird is going on with the volume, probably related to the above problem. Checkout a build from a few days ago and compare histograms. 
 * [x] Grid: Check consistency of centering when calculating histograms (appears to change result)
 * [x] General: Change system exits to throws

# Stuff to consider
## Grid:
 * Consider simply calculating and using the bounding box at initialization as the entire grid instead of wasting memory
 * Consider changing how it stores its atom data
 * Consider how to improve culling method
 * Consider removing all bounds checks
 * Consider using vectors instead of maps for storage. Currently insertion scales horribly with the number of atoms. 

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
 * ROOT (compile options: `cmake -DCMAKE_INSTALL_PREFIX=<install> CMAKE_CXX_STANDARD=17 -Dbuiltin_gsl=ON <source>`)
 * Boost
 * Elements

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
