/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/HydrationFitter.h>
#include <fitter/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <hist/Histogram.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <utility/Exceptions.h>
#include <settings/FitSettings.h>
#include <settings/HistogramSettings.h>
#include <io/ExistingFile.h>
#include <dataset/Dataset2D.h>
#include <mini/All.h>

using namespace fitter;

HydrationFitter::HydrationFitter(const SimpleDataset& data) : LinearFitter(data) {}
HydrationFitter::HydrationFitter(HydrationFitter&& other) : LinearFitter(std::move(other)), guess(std::move(other.guess)) {}
HydrationFitter::HydrationFitter(const io::ExistingFile& input) : LinearFitter(input) {}
HydrationFitter::HydrationFitter(const io::ExistingFile& input, std::unique_ptr<hist::ICompositeDistanceHistogram> h) : LinearFitter(input, std::move(h)) {initialize_guess();}
HydrationFitter::HydrationFitter(const SimpleDataset& data, std::unique_ptr<hist::ICompositeDistanceHistogram> h) : LinearFitter(data, std::move(h)) {initialize_guess();}
HydrationFitter::HydrationFitter(std::unique_ptr<hist::ICompositeDistanceHistogram> model) : HydrationFitter(std::move(model), Limit(settings::axes::qmin, settings::axes::qmax)) {}
HydrationFitter::HydrationFitter(std::unique_ptr<hist::ICompositeDistanceHistogram> model, const Limit& limits) : LinearFitter(std::move(model), limits) {initialize_guess();}

void HydrationFitter::initialize_guess() {
    auto lim = cast_h()->get_water_scaling_factor_limits();
    guess = mini::Parameter{"c", lim.center(), lim};
}

void HydrationFitter::set_algorithm(const mini::type& t) {fit_type = t;}
std::shared_ptr<FitResult> HydrationFitter::fit(const mini::type& algorithm) {
    fit_type = algorithm;
    return fit();
}

observer_ptr<hist::ICompositeDistanceHistogram> HydrationFitter::cast_h() const {
    return static_cast<hist::ICompositeDistanceHistogram*>(h.get());
}

std::shared_ptr<FitResult> HydrationFitter::fit() {
    mini::Parameter guess = {"c", 1, cast_h()->get_water_scaling_factor_limits()};
    std::function<double(std::vector<double>)> f = std::bind(&HydrationFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, std::move(f), guess, settings::fit::max_iterations);
    auto res = mini->minimize();

    // apply c
    cast_h()->apply_water_scaling_factor(res.get_parameter("c").value);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    std::shared_ptr<FitResult> ab_fit = fitter.fit();

    // update fitter object
    fitted = std::make_shared<FitResult>(res, res.fval, data.size()-1); // start with the fit performed here
    fitted->add_fit(ab_fit.get());                                      // add the a,b inner fit
    fitted->add_plots(this);                                            // make the result plottable
    fitted->evaluated_points = mini->get_evaluated_points();            // add the evaluated points
    return fitted;
}

double HydrationFitter::fit_chi2_only() {
    mini::Parameter guess = {"c", 1, cast_h()->get_water_scaling_factor_limits()};
    std::function<double(std::vector<double>)> f = std::bind(&HydrationFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, std::move(f), guess, settings::fit::max_iterations);
    auto res = mini->minimize();

    // apply c
    cast_h()->apply_water_scaling_factor(res.get_parameter("c").value);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    return fitter.fit_chi2_only();
}

FitResult::FitInfo HydrationFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("HydrationFitter::plot: Cannot plot before a fit has been made!");}

    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    cast_h()->apply_water_scaling_factor(c);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // calculate the scaled I model values
    std::vector<double> I_scaled(data.size()); // spliced data
    std::vector<double> ym_scaled(ym.size()); // original scaled data
    std::transform(Im.begin(), Im.end(), I_scaled.begin(), [&a, &b] (double I) {return I*a+b;});
    std::transform(ym.begin(), ym.end(), ym_scaled.begin(), [&a, &b] (double I) {return I*a+b;});

    FitResult::FitInfo graphs;
    graphs.fitted_intensity_interpolated = SimpleDataset(data.x(), I_scaled);
    graphs.fitted_intensity = SimpleDataset(h->get_q_axis(), ym_scaled);
    graphs.dataset = SimpleDataset(data.x(), data.y(), data.yerr());

    auto lim = graphs.dataset.get_xlimits();
    lim.expand(0.05);
    graphs.fitted_intensity.limit_x(lim);
    return graphs;
}

SimpleDataset HydrationFitter::plot_residuals() {
    if (fitted == nullptr) {throw except::bad_order("HydrationFitter::plot_residuals: Cannot plot before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    cast_h()->apply_water_scaling_factor(c);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // calculate the residuals
    std::vector<double> residuals(data.size());
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        residuals[i] = ((data.y(i) - (a*Im[i]+b))/data.yerr(i));
    }

    std::vector<double> xerr(data.size(), 0);
    return Dataset2D(data.x(), std::move(residuals), std::move(xerr), data.yerr());
}

double HydrationFitter::chi2(const std::vector<double>& params) {
    double c = params[0];

    // apply c
    cast_h()->apply_water_scaling_factor(c);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    auto[a, b] = fitter.fit_params_only();

    // calculate chi2
    double chi = 0;
    for (int i = 0; i < static_cast<int>(data.size()); i++) {
        double v = (data.y(i) - (a*Im[i]+b))/data.yerr(i);
        chi += v*v;
    }

    return chi;
}

unsigned int HydrationFitter::dof() const {
    return data.size() - 4;
}

double HydrationFitter::get_intercept() {
    if (fitted == nullptr) {throw except::bad_order("HydrationFitter::get_intercept: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    cast_h()->apply_water_scaling_factor(c);
    std::vector<double> ym = h->debye_transform().get_counts();
    math::CubicSpline s(h->get_q_axis(), ym);
    return a*s.spline(0) + b;
}

SimpleDataset HydrationFitter::get_model_dataset() {
    if (fitted == nullptr) {throw except::bad_order("HydrationFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    cast_h()->apply_water_scaling_factor(c);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);
    std::transform(Im.begin(), Im.end(), Im.begin(), [&a, &b] (double I) {return I*a+b;});

    return SimpleDataset(data.x(), Im, "q", "I"); 
}

SimpleDataset HydrationFitter::get_model_dataset(const std::vector<double>& q) {
    if (fitted == nullptr) {throw except::bad_order("HydrationFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
 
    auto model = get_model_dataset();
    return model.interpolate(q);
}

SimpleDataset HydrationFitter::get_dataset() const {
    return data;
}

void HydrationFitter::set_scattering_hist(std::unique_ptr<hist::ICompositeDistanceHistogram> h) {
    this->h = std::move(h);
}

observer_ptr<hist::ICompositeDistanceHistogram> HydrationFitter::get_scattering_hist() {
    return static_cast<hist::ICompositeDistanceHistogram*>(h.get());
}

void HydrationFitter::set_guess(const mini::Parameter& guess) {
    this->guess = guess;
}

void HydrationFitter::operator=(HydrationFitter&& other) {
    LinearFitter::operator=(std::move(other));
    guess = std::move(other.guess);
}