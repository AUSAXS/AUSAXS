# Todo
 * [ ] SLICES: change *Matrix pointers to &vector<double> for better efficiency
 * [ ] ATOM: Const uid 
 * [ ] CONSTANTS: Rethink how to determine charge densities for arbitrary ligands

# Stuff to consider
## Grid:
 * Consider simply calculating and using the bounding box at initialization as the entire grid instead of wasting memory
 * Consider changing how it stores its atom data
 * Consider how to improve culling method

## ScatteringHistogram:
 * OPTIMIZATION: Optional argument of q-values to calculate I(q) for - this would remove the necessity of splicing in the IntensityFitter

# Dependencies
Maybe bundle them somehow to make it easier to install?
 * ROOT
 * Boost

# Other personal notes
## Articles
 * B. Jacrot rep prog phys ~ 1980
 * Peter Zipper ~1980
 * ATSAS Crysol
