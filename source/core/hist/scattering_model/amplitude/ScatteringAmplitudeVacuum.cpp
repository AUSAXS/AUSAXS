// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/amplitude/ScatteringAmplitudeVacuum.h>

#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/ScatteringModelNoExv.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::hist;

ScatteringAmplitudeVacuum::ScatteringAmplitudeVacuum(std::vector<hist::detail::CompactCoordinates>&& groups) : ScatteringAmplitude(std::move(groups)) {}

container::Container2D<double> ScatteringAmplitudeVacuum::get_excess_form_factors(const std::vector<double>& q, const ScatteringParameters& params) const {
    const auto& ff_indices = form_factor::manager::get_active_indices();
    const int n_ff = form_factor::get_active_count();
    const int n_q = static_cast<int>(q.size());

    // the atoms are damped by the atomic Debye-Waller factor, while the waters are scaled by cw
    std::vector<double> atomic_scale(n_q), water_scale(n_q, params.cw);
    std::ranges::transform(q, atomic_scale.begin(), [&params] (double qk) {return ScatteringModelNoExv::get_atomic_debye_waller_factor(qk, params.DW_sigma_atomic);});

    container::Container2D<double> f(n_ff, n_q, 0); // no scatterer has the excluded volume slot
    for (int ff = form_factor::start_index_for_explicit_exv(); ff < n_ff; ++ff) {
        auto type = static_cast<form_factor::form_factor_t>(ff_indices[ff]);
        const auto& scale = ff == form_factor::water_bin ? water_scale : atomic_scale;
        std::ranges::transform(q, scale, f.begin(ff), [type] (double qk, double s) {return s*form_factor::manager::evaluate(type, qk);});
    }
    return f;
}
