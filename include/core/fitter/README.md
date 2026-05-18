This folder contains the fitters that match a calculated scattering profile against experimental data.

Contents
- `Fitter.h` — the common fitter interface.
- `LinearFitter.h` — fits the linear parameters (scaling and background) by least squares.
- `SmartFitter.h` — additionally optimizes the non-linear parameters such as excluded volume and hydration-shell scaling.
- `FitResult.h`, `FitReporter.h` — store the outcome of a fit and present it to the user.
