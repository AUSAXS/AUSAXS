// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <fitter/detail/LinearLeastSquares.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <math/CubicSpline.h>

#include <cmath>
#include <cassert>

using namespace ausaxs::fitter::detail;

LinearLeastSquares::LinearLeastSquares(const std::vector<double>& data, const std::vector<double>& model) : data(data), model(model), inv_sigma(this->data.size(), 1) {
    assert(this->data.size() == this->model.size() && "LinearFitter::LinearFitter: Data and model must have the same size.");
}

LinearLeastSquares::LinearLeastSquares(const std::vector<double>& data, const std::vector<double>& model, const std::vector<double>& errors) 
    : LinearLeastSquares(data, model)
{
    assert(data.size() == errors.size() && "LinearFitter::LinearFitter: Data and errors must have the same size.");
    for (int i = 0; i < static_cast<int>(errors.size()); ++i) {
        inv_sigma[i] = 1./errors[i];
    }
}

std::vector<double> LinearLeastSquares::fit_params_only() {
    double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        double inv_sig2 = inv_sigma[i]*inv_sigma[i];
        S += inv_sig2;
        Sx += data[i]*inv_sig2;
        Sy += model[i]*inv_sig2;
        Sxx += std::pow(data[i], 2)*inv_sig2;
        Sxy += data[i]*model[i]*inv_sig2;
    }

    double delta = S*Sxx - Sx*Sx;
    assert(delta != 0 && "LinearLeastSquares::fit_params_only: Division by zero.");
    double a = (S*Sxy - Sx*Sy)/delta;
    double b = (Sxx*Sy - Sx*Sxy)/delta;
    double a_err = S/delta;
    double b_err = Sxx/delta;
    return {a, b, a_err, b_err};
}

std::unique_ptr<ausaxs::fitter::FitResult> LinearLeastSquares::fit() {
    auto p = fit_params_only();

    std::unique_ptr<FitResult> f = std::make_unique<FitResult>();
    f->parameters = {{"a", p[0], std::sqrt(p[2])}, {"b", p[1], std::sqrt(p[3])}};
    f->dof = dof();
    f->fval = chi2(p);
    f->fevals = 1;
    return f;
}

std::vector<double> LinearLeastSquares::get_model_curve(const std::vector<double>& p) {
    assert(p.size() == 2 && "LinearFitter::get_model_curve: Invalid number of parameters.");
    std::vector<double> model(data.size());
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        model[i] = p[0]*data[i] + p[1];
    }
    return model;
}

std::vector<double> LinearLeastSquares::get_model_curve() {
    auto p = fit_params_only();
    return get_model_curve({p[0], p[1]});
}

std::vector<double> LinearLeastSquares::get_residuals(const std::vector<double>& p) {
    assert(2 <= p.size() && "LinearFitter::get_residuals: Invalid number of parameters. Expected at least 2 (a, b).");
    std::vector<double> residuals(data.size());
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        residuals[i] = (model[i] - (p[0]*data[i] + p[1]))*inv_sigma[i];
    }
    return residuals;
}

unsigned int LinearLeastSquares::dof() const {return data.size() - 2;}

unsigned int LinearLeastSquares::size() const {return data.size();}