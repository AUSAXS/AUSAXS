#include <fitter/refactored/Fitter.h>

#include <numeric>

double fitter::Fitter::chi2(const std::vector<double>& params) {
    auto residuals = get_residuals(params);
    return std::accumulate(residuals.begin(), residuals.end(), 0.0, [] (double sum, double val) {return sum + val*val;});
}