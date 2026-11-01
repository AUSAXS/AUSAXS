// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/IScatteringModelExv.h>

using namespace ausaxs;

Limit hist::IScatteringModelExv::get_excluded_volume_scaling_factor_limits() const {return {0.92, 1.08};}

Limit hist::IScatteringModelExv::get_solvent_density_scaling_factor_limits() const {return {0.5, 2};}

Limit hist::IScatteringModelExv::get_debye_waller_factor_limits() const {return {0.0, 5};}

std::vector<double>& hist::IScatteringModelExv::get_total_raw_counts() {
    return const_cast<std::vector<double>&>(const_cast<const IScatteringModelExv*>(this)->get_total_raw_counts());
}

const std::vector<double>& hist::IScatteringModelExv::get_raw_counts() const {
    return get_total_raw_counts();
}

std::vector<double>& hist::IScatteringModelExv::get_raw_counts() {
    return get_total_raw_counts();
}