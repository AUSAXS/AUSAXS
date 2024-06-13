/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/SimpleLeastSquares.h>
#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <fitter/FitPlots.h>
#include <utility/Exceptions.h>

#include <cmath>

using namespace fitter;

SimpleLeastSquares::SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model) : data(data, model, std::vector<double>(data.size(), 1)) {}

SimpleLeastSquares::SimpleLeastSquares(const std::vector<double> data, const std::vector<double> model, const std::vector<double> errors) : data(data, model, errors) {}

SimpleLeastSquares::SimpleLeastSquares(const SimpleDataset& data) : data(data) {}

SimpleLeastSquares::SimpleLeastSquares(SimpleDataset&& data) : data(std::move(data)) {}

std::pair<double, double> SimpleLeastSquares::fit_params_only() {
    S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for (unsigned i = 0; i < data.size(); ++i) {
        double sig2 = pow(data.yerr(i), 2);
        S += 1./sig2;
        Sx += data.x(i)/sig2;
        Sy += data.y(i)/sig2;
        Sxx += pow(data.x(i), 2)/sig2;
        Sxy += data.x(i)*data.y(i)/sig2;
    }

    delta = S*Sxx - Sx*Sx;
    a = (S*Sxy - Sx*Sy)/delta;
    b = (Sxx*Sy - Sx*Sxy)/delta;
    return std::make_pair(a, b);
}

double SimpleLeastSquares::fit_chi2_only() {
    [[maybe_unused]] auto res = fit_params_only();
    return chi2({});
}

std::shared_ptr<Fit> SimpleLeastSquares::fit() {
    if (delta == 0) {[[maybe_unused]] auto res = fit_params_only();}
    double a_err2 = S/delta; // squared sigmas
    double b_err2 = Sxx/delta; 

    // double cov_ab = -Sx/delta;
    // double Q = ROOT::Math::inc_gamma((double) data.size()/2 -1, chi2()/2);

    std::shared_ptr<Fit> f = std::make_shared<Fit>();
    f->parameters = {{"a", a, sqrt(a_err2)}, {"b", b, sqrt(b_err2)}};
    f->dof = data.size() - 2;
    f->fval = chi2({});
    f->fevals = 1;
    // f->status = Q > 0.001 ? 0 : 1;
    return f;
}

double SimpleLeastSquares::chi2(const std::vector<double>&) {
    double chi = 0;
    for (unsigned int i = 0; i < data.size(); ++i) {
        chi += pow((data.y(i) - (a*data.x(i) + b))/data.yerr(i), 2);
    }

    return chi;
}

FitPlots SimpleLeastSquares::plot() {
    throw except::unexpected("SimpleLeastSquares::plot: Not implemented yet. ");
}

SimpleDataset SimpleLeastSquares::plot_residuals() {
    throw except::unexpected("SimpleLeastSquares::plot_residuals: Not implemented yet. ");
}

unsigned int SimpleLeastSquares::dof() const {return data.size() - 2;}

unsigned int SimpleLeastSquares::size() const {return data.size();}

std::shared_ptr<Fit> SimpleLeastSquares::get_fit() const {
    throw except::unexpected("SimpleLeastSquares::get_fit: Not implemented yet. ");
}