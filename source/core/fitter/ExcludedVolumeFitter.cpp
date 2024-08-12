/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <hist/Histogram.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <utility/Exceptions.h>
#include <mini/All.h>
#include <settings/GeneralSettings.h>

using namespace fitter;

ExcludedVolumeFitter::ExcludedVolumeFitter(const SimpleDataset& data) : HydrationFitter(data) {}

ExcludedVolumeFitter::ExcludedVolumeFitter(const SimpleDataset& saxs, std::unique_ptr<hist::ICompositeDistanceHistogram> h) : ExcludedVolumeFitter(saxs) {
    set_scattering_hist(std::move(h));
}

void ExcludedVolumeFitter::validate_histogram() const {
    if (h == nullptr) {throw except::invalid_argument("ExcludedVolumeFitter::validate_histogram: Cannot fit without a histogram.");}
    if (dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h.get()) == nullptr) {
        throw except::invalid_argument("ExcludedVolumeFitter::validate_histogram: Histogram must support excluded volume operations.");
    }
}

void ExcludedVolumeFitter::initialize_guess() {
    auto lim = cast_exv()->get_excluded_volume_scaling_factor_limits();
    guess_exv = mini::Parameter{"d", lim.min < 1 && 1 < lim.max ? 1 : lim.center(), lim};
}

std::shared_ptr<FitResult> ExcludedVolumeFitter::fit() {
    fit_type = mini::type::DEFAULT;
    settings::general::verbose = false;
    std::function<double(std::vector<double>)> f = std::bind(&ExcludedVolumeFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, std::move(f), {guess, guess_exv});
    auto res = mini->minimize();

    update_excluded_volume(res.get_parameter("d").value);
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

    settings::general::verbose = true;
    return fitted;
}

double ExcludedVolumeFitter::fit_chi2_only() {
    fit_type = mini::type::DEFAULT;
    settings::general::verbose = false;
    std::function<double(std::vector<double>)> f = std::bind(&ExcludedVolumeFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, std::move(f), {guess, guess_exv});
    auto res = mini->minimize();

    update_excluded_volume(res.get_parameter("d").value);
    cast_h()->apply_water_scaling_factor(res.get_parameter("c").value);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    return fitter.fit_chi2_only();
}

FitResult::FitInfo ExcludedVolumeFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("ExcludedVolumeFitter::plot: Cannot plot before a fit has been made!");}
    update_excluded_volume(fitted->get_parameter("d").value);
    return HydrationFitter::plot();
}

SimpleDataset ExcludedVolumeFitter::plot_residuals() {
    if (fitted == nullptr) {throw except::bad_order("ExcludedVolumeFitter::plot_residuals: Cannot plot before a fit has been made!");}
    update_excluded_volume(fitted->get_parameter("d").value);
    return HydrationFitter::plot_residuals();
}

double ExcludedVolumeFitter::chi2(const std::vector<double>& params) {
    update_excluded_volume(params[1]);
    return HydrationFitter::chi2(params);
}

double ExcludedVolumeFitter::get_intercept() {
    if (fitted == nullptr) {throw except::bad_order("HydrationFitter::get_intercept: Cannot determine model intercept before a fit has been made!");}
    update_excluded_volume(fitted->get_parameter("d").value);
    return HydrationFitter::get_intercept();
}

void ExcludedVolumeFitter::update_excluded_volume(double d) {
    #if defined(DEBUG)
        auto res = dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h.get());
        if (res == nullptr) {throw except::invalid_operation("ExcludedVolumeFitter::update_excluded_volume: Cannot cast histogram to ICompositeDistanceHistogramExv!");}
    #endif
    static_cast<hist::ICompositeDistanceHistogramExv*>(h.get())->apply_excluded_volume_scaling_factor(d);
}

SimpleDataset ExcludedVolumeFitter::get_model_dataset() {
    if (fitted == nullptr) {throw except::bad_order("ExcludedVolumeFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
    update_excluded_volume(fitted->get_parameter("d").value);
    return HydrationFitter::get_model_dataset();
}

SimpleDataset ExcludedVolumeFitter::get_model_dataset(const std::vector<double>& q) {
    if (fitted == nullptr) {throw except::bad_order("ExcludedVolumeFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
    update_excluded_volume(fitted->get_parameter("d").value);
    return HydrationFitter::get_model_dataset(q);
}

SimpleDataset ExcludedVolumeFitter::get_dataset() const {
    return data;
}

unsigned int ExcludedVolumeFitter::dof() const {
    return data.size() - 5;
}

void ExcludedVolumeFitter::set_guess(mini::Parameter guess_hydration, mini::Parameter guess_exv) {
    this->guess = std::move(guess_hydration);
    this->guess_exv = std::move(guess_exv);
}

observer_ptr<hist::ICompositeDistanceHistogramExv> ExcludedVolumeFitter::cast_exv() const {
    return static_cast<hist::ICompositeDistanceHistogramExv*>(h.get());
}

void ExcludedVolumeFitter::set_scattering_hist(std::unique_ptr<hist::ICompositeDistanceHistogram> h) {
    HydrationFitter hfit(data, std::move(h));
    auto hres = hfit.fit();
    double c = hres->get_parameter("c").value;
    HydrationFitter::operator=(std::move(hfit));
    initialize_guess();
    validate_histogram();
    guess.guess = c;
}