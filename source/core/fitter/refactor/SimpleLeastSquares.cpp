/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/refactored/SimpleLeastSquares.h>
#include <fitter/FitResult.h>

#include <cmath>
#include <cassert>

using namespace fitter;

SimpleLeastSquares::SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model) : data(data), model(model), inv_sigma2(data.size(), 1) {
    assert(data.size() == model.size() && "SimpleLeastSquares::SimpleLeastSquares: Data and model must have the same size.");
}

SimpleLeastSquares::SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model, const std::vector<double> errors) : SimpleLeastSquares(data, model) {
    assert(data.size() == errors.size() && "SimpleLeastSquares::SimpleLeastSquares: Data and errors must have the same size.");
    for (unsigned i = 0; i < errors.size(); ++i) {
        inv_sigma2[i] = 1./pow(errors[i], 2);
    }
}

std::array<double, 4> SimpleLeastSquares::fit_params_only() {
    double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for (unsigned i = 0; i < data.size(); ++i) {
        S += inv_sigma2[i];
        Sx += data[i]*inv_sigma2[i];
        Sy += model[i]*inv_sigma2[i];
        Sxx += std::pow(data[i], 2)*inv_sigma2[i];
        Sxy += data[i]*model[i]*inv_sigma2[i];
    }

    double delta = S*Sxx - Sx*Sx;
    double a = (S*Sxy - Sx*Sy)/delta;
    double b = (Sxx*Sy - Sx*Sxy)/delta;
    double a_err = S/delta;
    double b_err = Sxx/delta;
    return {a, b, a_err, b_err};
}

double SimpleLeastSquares::fit_chi2_only() {
    std::ignore = fit_params_only(); // update a and b
    return chi2({});
}

std::unique_ptr<FitResult> SimpleLeastSquares::fit() {
    auto[a, b, aerr2, berr2] = fit_params_only();

    std::unique_ptr<FitResult> f = std::make_unique<FitResult>();
    f->parameters = {{"a", a, std::sqrt(aerr2)}, {"b", b, std::sqrt(berr2)}};
    f->dof = data.size() - 3;
    f->fval = chi2({});
    f->fevals = 1;
    return f;
}

double SimpleLeastSquares::chi2(const std::vector<double>&) {
    auto[a, b, _, __] = fit_params_only();
    double chi = 0;
    for (unsigned int i = 0; i < data.size(); ++i) {
        chi += std::pow((model[i] - (a*data[i] + b)), 2)*inv_sigma2[i];
    }

    return chi;
}

unsigned int SimpleLeastSquares::dof() const {return data.size() - 2;}

unsigned int SimpleLeastSquares::size() const {return data.size();}