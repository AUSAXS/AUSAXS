This folder contains the function minimizers used for parameter optimization, for example during fitting.

A range of strategies is available — golden-section search, plain and limited scans, a minimum explorer, and a wrapper around the dlib optimizers. `MinimizerFactory.h` selects an implementation by name, and `All.h` is a convenience header pulling in every minimizer.
