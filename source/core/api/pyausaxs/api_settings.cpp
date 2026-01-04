// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_settings.h>
#include <settings/SettingsIO.h>
#include <settings/All.h>

using namespace ausaxs;

void set_exv_settings(
    const char* exv_model,
    int* status
) {return execute_with_catch([&]() {
    settings::detail::parse_option("exv_model", {std::string(exv_model)});
}, status);}

void set_fit_settings(
    unsigned int N,
    unsigned int max_iterations,
    bool fit_excluded_volume,
    bool fit_solvent_density,
    bool fit_hydration,
    bool fit_atomic_debye_waller,
    bool fit_exv_debye_waller,
    int* status
) {return execute_with_catch([&]() {
    settings::fit::N = N;
    settings::fit::max_iterations = max_iterations;
    settings::fit::fit_excluded_volume = static_cast<bool>(fit_excluded_volume);
    settings::fit::fit_solvent_density = static_cast<bool>(fit_solvent_density);
    settings::fit::fit_hydration = static_cast<bool>(fit_hydration);
    settings::fit::fit_atomic_debye_waller = static_cast<bool>(fit_atomic_debye_waller);
    settings::fit::fit_exv_debye_waller = static_cast<bool>(fit_exv_debye_waller);
}, status);}

void set_grid_settings(
    double water_scaling,
    double cell_width,
    double scaling,
    double min_exv_radius,
    unsigned int min_bins,
    int* status
) {return execute_with_catch([&]() {
    settings::grid::water_scaling = water_scaling;
    settings::grid::cell_width = cell_width;
    settings::grid::scaling = scaling;
    settings::grid::min_exv_radius = min_exv_radius;
    settings::grid::min_bins = min_bins;
}, status);}

void set_hist_settings(
    unsigned int skip,
    double qmin,
    double qmax,
    bool weighted_bins,
    double bin_width,
    unsigned int bin_count,
    int* status
) {return execute_with_catch([&]() {
    settings::axes::skip = skip;
    settings::axes::qmin = qmin;
    settings::axes::qmax = qmax;
    settings::axes::bin_width = bin_width;
    settings::axes::bin_count = bin_count;
    settings::hist::weighted_bins = weighted_bins;
}, status);}

void set_molecule_settings(
    bool center,
    bool throw_on_unknown_atom,
    bool implicit_hydrogens,
    bool use_occupancy,
    const char* exv_set,
    const char* hydration_strategy,
    int* status
) {return execute_with_catch([&]() {
    settings::molecule::center = center;
    settings::molecule::throw_on_unknown_atom = throw_on_unknown_atom;
    settings::molecule::implicit_hydrogens = implicit_hydrogens;
    settings::molecule::use_occupancy = use_occupancy;
    settings::detail::parse_option("exv_volume", {std::string(exv_set)});
    settings::detail::parse_option("hydration_strategy", {std::string(hydration_strategy)});
}, status);}

void set_general_settings(
    bool offline,
    bool verbose,
    bool warnings,
    unsigned int threads,
    int* status
) {return execute_with_catch([&]() {
    settings::general::offline = offline;
    settings::general::verbose = verbose;
    settings::general::warnings = warnings;
    settings::general::threads = threads;
}, status);}