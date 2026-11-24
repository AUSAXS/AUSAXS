// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/amplitude/ScatteringAmplitudeVacuum.h>

#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/ScatteringModelNoExv.h>

using namespace ausaxs;
using namespace ausaxs::hist;

ScatteringAmplitudeVacuum::ScatteringAmplitudeVacuum(std::vector<hist::detail::CompactCoordinates>&& groups) : ScatteringAmplitude(std::move(groups)) {}

std::vector<double> ScatteringAmplitudeVacuum::get_excess_form_factors(double q, const ScatteringParameters& params) const {
    const auto& ff_indices = form_factor::manager::get_active_indices();
    const int n_ff = form_factor::get_active_count();
    const double B = ScatteringModelNoExv::get_atomic_debye_waller_factor(q, params.DW_sigma_atomic);
    constexpr int water = form_factor::water_bin;

    std::vector<double> f(n_ff, 0); // no scatterer has the excluded volume slot
    for (int ff = form_factor::start_index_for_explicit_exv(); ff < n_ff; ++ff) {
        if (ff == water) {continue;}
        f[ff] = B*form_factor::manager::evaluate(static_cast<form_factor::form_factor_t>(ff_indices[ff]), q);
    }
    f[water] = params.cw*form_factor::manager::evaluate(static_cast<form_factor::form_factor_t>(ff_indices[water]), q);
    return f;
}
