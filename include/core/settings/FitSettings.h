#pragma once

namespace ausaxs::settings::fit {
    extern bool verbose;                 // Decides if the fitting process will be verbose.
    extern unsigned int N;               // Number of points sampled when discretizing a model scattering curve
    extern unsigned int max_iterations;  // Maximum number of iterations in the fitting process
    extern bool fit_excluded_volume;     // Enable fitting of the excluded volume solvent.
    extern bool fit_solvent_density;     // Enable fitting of the solvent density for the excluded volume.
    extern bool fit_hydration;           // Enable fitting of the hydration shell.
}