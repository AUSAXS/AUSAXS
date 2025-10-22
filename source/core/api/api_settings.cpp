#include <api/api_settings.h>

void set_exv_settings(
    const char* exv_model,
    int* status
);

void set_fit_settings(
    unsigned int N,
    unsigned int max_iterations,
    int fit_excluded_volume,
    int fit_solvent_density,
    int fit_hydration,
    int fit_atomic_debye_waller,
    int fit_exv_debye_waller,
    int* status
);

void set_grid_settings(
    double water_scaling,
    double cell_width,
    double scaling,
    double min_exv_radius,
    unsigned int min_bins,
    int* status
);

void set_hist_settings(
    unsigned int skip,
    double qmin,
    double qmax,
    bool weighted_bins
);

void set_molecule_settings(
    bool center,
    bool throw_on_unknown_atom,
    bool implicit_hydrogens,
    bool use_occupancy,
    const char* exv_set,
    const char* hydration_strategy
);