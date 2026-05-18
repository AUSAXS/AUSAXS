This folder contains the standalone math library used across all of AUSAXS.

It provides the basic linear-algebra types `Vector`, `Vector3`, and `Matrix`, a range of matrix decompositions and solvers (LUP, QR, Givens), and assorted numerical tools such as cubic-spline interpolation, a peak finder, a moving averager, and statistics helpers. `indexers/` and `slices/` offer views into multidimensional data. The library has no dependencies on the rest of the codebase.
