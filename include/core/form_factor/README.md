This folder contains the atomic form factors that describe how each atom type scatters X-rays as a function of the scattering vector q.

Contents
- `FormFactor.h`, `FormFactorType.h` — the form-factor class and the enumeration of supported types.
- `FormFactorTable.h` — tabulated form-factor coefficients.
- `ExvFormFactor.h`, `ExvTable.h` — excluded-volume form factors, representing the solvent displaced by each atom.
- `NormalizedFormFactor.h` — form factors normalized for use in the histogram-based intensity calculation.
- `FormFactorConcepts.h` — C++ concepts constraining form-factor template parameters.
- `lookup/` — the manager and product tables that cache form-factor products for every pair of types, avoiding repeated evaluation during histogram-to-intensity conversion.
