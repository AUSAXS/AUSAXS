This folder contains the function minimizers used for parameter optimization, for example during fitting.

Every minimizer takes the function to minimize as a vector of residuals r(p), and minimizes chi2 = sum_i r_i(p)^2. Levenberg-Marquardt (the default) works on the residuals directly, while the others only use their sum.

A range of strategies is available — Levenberg-Marquardt, golden-section search, plain and limited scans, a minimum explorer, and a wrapper around the dlib optimizers. `MinimizerFactory.h` selects an implementation by name, and `All.h` is a convenience header pulling in every minimizer.
