/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/Fitter.h>

#include <numeric>

using namespace ausaxs;

double fitter::Fitter::chi2(const std::vector<double>& params) {
    auto residuals = get_residuals(params);
    double chi2 = std::accumulate(residuals.begin(), residuals.end(), 0.0, [] (double sum, double val) {return sum + val*val;});
    std::cout << "Chi2: " << chi2 << std::endl;
    return std::accumulate(residuals.begin(), residuals.end(), 0.0, [] (double sum, double val) {return sum + val*val;});
}

double fitter::Fitter::fit_chi2_only() {
    return chi2(fit_params_only());
}