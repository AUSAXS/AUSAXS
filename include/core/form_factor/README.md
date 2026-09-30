This folder contains the atomic form factors that describe how each atom type scatters X-rays or neutrons as a function of the scattering vector q. The probe is selected at runtime by `settings::scattering::radiation`, and the probe-specific code lives in the `form_factor::xray` and `form_factor::neutron` namespaces.

Contents
- `FormFactorType.h` — the enumeration of supported types, and `ff_info_table`, the descriptor table holding each type's name, element, hydrogens, electrons, mass and X-ray coefficients. To add a type, add it to the enum and give it a row there.
- `FormFactor.h` — the X-ray form-factor class `xray::FormFactor`, and the lookup `xray::raw` built from the descriptor table.
- `FormFactorTable.h` — tabulated five-Gaussian X-ray form-factor coefficients, with their literature sources.
- `NeutronFormFactor.h` — the neutron form-factor class `neutron::FormFactor`, holding both the orientationally averaged amplitude of a group and its self-term, and the lookups `neutron::protonated` and `neutron::deuterated`. The table itself is generated at compile time in `NeutronFormFactor.cpp` from the coherent scattering lengths and the X-H geometry.
- `ExvTable.h` — `exv_info_table`, the optional excluded-volume descriptor table with one row per type and one column per displaced-volume set. A type missing from the current set cannot be used with the Fraser-based models and is treated as OTHER.
- `ExvFormFactor.h` — excluded-volume form factors, representing the solvent displaced by each atom.
- `NormalizedFormFactor.h` — X-ray form factors normalized to 1 at q = 0.
- `FormFactorConcepts.h` — C++ concepts constraining form-factor template parameters.
- `lookup/` — the manager and product tables that cache form-factor products for every pair of types, avoiding repeated evaluation during histogram-to-intensity conversion.
