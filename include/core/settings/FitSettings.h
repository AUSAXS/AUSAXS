#pragma once

namespace ausaxs::settings {
    struct fit {
        static bool verbose;                 // Decides if the fitting process will be verbose.
        static unsigned int N;               // Number of points sampled when discretizing a model scattering curve
        static unsigned int max_iterations;  // Maximum number of iterations in the fitting process
        static bool fit_excluded_volume;     // Enable fitting of the excluded volume solvent.
        static bool fit_solvent_density;     // Enable fitting of the solvent density for the excluded volume.
        static bool fit_hydration;           // Enable fitting of the hydration shell.
        static bool fit_atomic_debye_waller; // Enable fitting of the atomic form factor debye-waller factor.
        static bool fit_exv_debye_waller;    // Enable fitting of the excluded volume form factor debye-waller factor.
    };
}