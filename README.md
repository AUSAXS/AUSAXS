# Todo
 * [ ] ATOM: Const uid 
 * [ ] CONSTANTS: Rethink how to determine charge densities for arbitrary ligands
 * [ ] GRID: Fix AxesPlacement - exact water locations, iterate over grid

# Stuff to consider
## Grid:
 * Consider simply calculating and using the bounding box at initialization as the entire grid instead of wasting memory
 * Consider changing how it stores its atom data
 * Consider how to improve culling method
 * Consider removing all bounds checks

## ScatteringHistogram:
 * Optional argument of q-values to calculate I(q) for - this would remove the necessity of splicing in the IntensityFitter
 * Take a closer look at the form factor
 * Convert `a` to a more sensible output (units)

## Body:
 * Consider removing hydration_atoms, they're not supposed to be used anyway

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
