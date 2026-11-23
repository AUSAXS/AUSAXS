// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringAmplitudeFraser.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/ScatteringModelFraser.h>

using namespace ausaxs;
using namespace ausaxs::hist;

std::vector<double> ScatteringAmplitudeFraser::form_factors(double q, const ScatteringParameters& params) const {
    auto f = ScatteringAmplitudeNoExv::form_factors(q, params);

    // the excluded volume tables of the intensity are only tabulated on the default q axis, so its form factors are evaluated directly
    const auto& ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
    const auto exv_set = form_factor::ExvTableManager::get_current_exv_form_factor_set();
    const double x = params.crho*ScatteringModelFraser::exv_factor(q, params.cx)
        *ScatteringModelFraser::get_exv_debye_waller_factor(q, params.DW_sigma_exv);
    for (int ff = form_factor::start_index_for_explicit_exv(); ff < form_factor::get_active_count(); ++ff) {
        if (ff == form_factor::water_bin) {continue;}
        auto type = static_cast<form_factor::form_factor_t>(ff_indices[ff]);
        if (!exv_set.contains(type)) {continue;} // a type without a volume has no excluded volume in the intensity either
        f[ff] -= x*exv_set.get(type).evaluate(q);
    }
    return f;
}
