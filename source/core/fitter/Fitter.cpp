// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <fitter/Fitter.h>

#include <functional>
#include <numeric>

using namespace ausaxs;

double fitter::Fitter::chi2(const std::vector<double>& params) {
    auto residuals = get_residuals(params);
    return std::transform_reduce(residuals.begin(), residuals.end(), 0.0, std::plus{}, [] (double val) {return val*val;});
}

double fitter::Fitter::fit_chi2_only() {
    return chi2(fit_params_only());
}