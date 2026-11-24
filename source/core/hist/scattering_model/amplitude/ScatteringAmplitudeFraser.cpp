// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/amplitude/ScatteringAmplitudeFraser.h>

#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/ScatteringModelFraser.h>

using namespace ausaxs;
using namespace ausaxs::hist;

container::Container2D<double> ScatteringAmplitudeFraser::get_excess_form_factors(const std::vector<double>& q, const ScatteringParameters& params) const {
    auto f = ScatteringAmplitudeVacuum::get_excess_form_factors(q, params);
    const int n_q = static_cast<int>(q.size());

    // the form factors are evaluated directly rather than read from the tables of the intensity, which keeps them exact at any q, including q = 0 and beyond the default q axis
    const auto& ff_indices = form_factor::manager::get_active_indices();
    const auto exv_set = form_factor::ExvTableManager::get_current_exv_form_factor_set();
    std::vector<double> exv_scale(n_q);
    for (int k = 0; k < n_q; ++k) {
        exv_scale[k] = params.crho*ScatteringModelFraser::exv_factor(q[k], params.cx)*ScatteringModelFraser::get_exv_debye_waller_factor(q[k], params.DW_sigma_exv);
    }

    for (int ff = form_factor::start_index_for_explicit_exv(); ff < form_factor::get_active_count(); ++ff) {
        if (ff == form_factor::water_bin) {continue;}
        auto type = static_cast<form_factor::form_factor_t>(ff_indices[ff]);
        if (!exv_set.contains(type)) {continue;} // a type without a volume has no excluded volume in the intensity either
        const auto& fx = exv_set.get(type);
        for (int k = 0; k < n_q; ++k) {f(ff, k) -= exv_scale[k]*fx.evaluate(q[k]);}
    }
    return f;
}
