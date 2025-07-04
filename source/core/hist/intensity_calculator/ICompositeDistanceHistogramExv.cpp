// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>

using namespace ausaxs;

Limit hist::ICompositeDistanceHistogramExv::get_excluded_volume_scaling_factor_limits() const {return {0.92, 1.08};}

Limit hist::ICompositeDistanceHistogramExv::get_solvent_density_scaling_factor_limits() const {return {0.5, 2};}

Limit hist::ICompositeDistanceHistogramExv::get_debye_waller_factor_limits() const {return {0.0, 5};}