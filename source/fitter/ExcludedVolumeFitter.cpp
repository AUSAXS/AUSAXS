#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/Fit.h>
#include <fitter/FitPlots.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <hist/Histogram.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <utility/Exceptions.h>
#include <mini/All.h>
#include <plots/All.h>
#include <settings/GeneralSettings.h>
#include <data/Molecule.h>

using namespace fitter;

ExcludedVolumeFitter::ExcludedVolumeFitter(const io::ExistingFile& input, std::unique_ptr<hist::ICompositeDistanceHistogram> h) : HydrationFitter(), fit_type(mini::type::BFGS) {
    HydrationFitter hfit(input, std::move(h));
    auto hres = hfit.fit();
    double c = hres->get_parameter("c").value;
    if (c == 0) {this->guess = {{"c", 0, {0, 1}},  {"d", 1, {0, 1.5}}};}
    else {this->guess = {{"c", c, {c*0.8, c*1.2}}, {"d", 1, {0, 1.5}}};}
    HydrationFitter::operator=(std::move(hfit));
}

std::shared_ptr<Fit> ExcludedVolumeFitter::fit() {
    fit_type = mini::type::DLIB_GLOBAL;
    settings::general::verbose = false;
    std::function<double(std::vector<double>)> f = std::bind(&ExcludedVolumeFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, f, guess);
    auto res = mini->minimize();

    update_excluded_volume(res.get_parameter("d").value);
    cast_h()->apply_water_scaling_factor(res.get_parameter("c").value);
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    std::shared_ptr<Fit> ab_fit = fitter.fit();

    // update fitter object
    fitted = std::make_shared<Fit>(res, res.fval, data.size()-2); // start with the fit performed here
    fitted->add_fit(ab_fit);                                      // add the a,b inner fit
    fitted->add_plots(*this);                                     // make the result plottable
    fitted->evaluated_points = mini->get_evaluated_points();      // add the evaluated points

    settings::general::verbose = true;
    return fitted;
}

double ExcludedVolumeFitter::fit_chi2_only() {
    fit_type = mini::type::DLIB_GLOBAL;
    settings::general::verbose = false;
    std::function<double(std::vector<double>)> f = std::bind(&ExcludedVolumeFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, f, guess);
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

FitPlots ExcludedVolumeFitter::plot() {
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

void ExcludedVolumeFitter::set_guess(std::vector<mini::Parameter> guess) {
    this->guess = guess;
}