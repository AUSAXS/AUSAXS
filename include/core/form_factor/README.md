This folder contains the atomic form factors that describe how each atom type scatters X-rays as a function of the scattering vector q.

Contents
- `FormFactorType.h` — the enumeration of supported types, and `ff_info_table`, the descriptor table holding each type's name, element, hydrogens, electrons, mass and coefficients. To add a type, add it to the enum and give it a row there.
- `FormFactor.h` — the form-factor class, and the lookup built from the descriptor table.
- `FormFactorTable.h` — tabulated five-Gaussian form-factor coefficients, with their literature sources.
- `ExvTable.h` — `exv_info_table`, the optional excluded-volume descriptor table with one row per type and one column per displaced-volume set. A type missing from the current set cannot be used with the Fraser-based models (Fraser, CRYSOL, Pepsi) and is treated as OTHER.
- `ExvFormFactor.h` — excluded-volume form factors, representing the solvent displaced by each atom.
- `NormalizedFormFactor.h` — form factors normalized for use in the histogram-based intensity calculation.
- `FormFactorConcepts.h` — C++ concepts constraining form-factor template parameters.
- `lookup/` — the manager and product tables that cache form-factor products for every pair of types, avoiding repeated evaluation during histogram-to-intensity conversion.
