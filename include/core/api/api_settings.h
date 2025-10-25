#pragma once

#include <api/api_helper.h>

extern "C" API void set_general_settings(
    bool offline,
    bool verbose,
    bool warnings,
    unsigned int threads,
    int* status
);

extern "C" API void set_exv_settings(
    const char* exv_model,
    int* status
);

extern "C" API void set_fit_settings(
    unsigned int N,
    unsigned int max_iterations,
    bool fit_excluded_volume,
    bool fit_solvent_density,
    bool fit_hydration,
    bool fit_atomic_debye_waller,
    bool fit_exv_debye_waller,
    int* status
);

extern "C" API void set_grid_settings(
    double water_scaling,
    double cell_width,
    double scaling,
    double min_exv_radius,
    unsigned int min_bins,
    int* status
);

extern "C" API void set_hist_settings(
    unsigned int skip,
    double qmin,
    double qmax,
    bool weighted_bins,
    int* status
);

extern "C" API void set_molecule_settings(
    bool center,
    bool throw_on_unknown_atom,
    bool implicit_hydrogens,
    bool use_occupancy,
    const char* exv_set,
    const char* hydration_strategy,
    int* status
);