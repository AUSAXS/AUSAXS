This folder contains everything related to distance histograms — the intermediate representation used to evaluate the Debye equation efficiently.

Rather than summing over every atom pair directly, AUSAXS first bins the pairwise distances into a histogram, then transforms that histogram into a scattering profile. This is the core approximation that makes the calculations fast.

Contents
- `Histogram.h`, `Histogram2D.h` — basic histogram containers.
- `distribution/` — the 1D/2D/3D distance distributions, including weighted variants that retain the mean distance within each bin for improved accuracy.
- `distance_calculator/` — routines that compute pairwise distances and accumulate them into distributions.
- `histogram_manager/` — managers orchestrating histogram construction, including multithreaded and partial (incremental) variants used during optimization.
- `intensity_calculator/` — converts distance histograms into scattering intensities, with composite variants handling excluded volume and per-form-factor contributions. Includes interfaces to external conventions (`crysol`, `foxs`, `pepsi`).
