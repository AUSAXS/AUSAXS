#include <fitter/Fitter.h>

#include <numeric>

using namespace ausaxs;

double fitter::Fitter::chi2(const std::vector<double>& params) {
    auto residuals = get_residuals(params);
    return std::accumulate(residuals.begin(), residuals.end(), 0.0, [] (double sum, double val) {return sum + val*val;});
}

double fitter::Fitter::fit_chi2_only() {
    return chi2(fit_params_only());
}