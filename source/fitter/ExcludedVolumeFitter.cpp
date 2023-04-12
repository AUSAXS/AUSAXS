#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/Fit.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <hist/ScatteringHistogram.h>
#include <utility/Exceptions.h>
#include <mini/all.h>
#include <plots/all.h>

using namespace fitter;

ExcludedVolumeFitter::ExcludedVolumeFitter(std::string input, Protein& protein) : HydrationFitter(input, protein.get_histogram()), protein(protein) {
    HydrationFitter hfit(input, protein.get_histogram());
    auto hres = hfit.fit();
    double c = hres->get_parameter("c").value;
    this->guess = {{"c", c, {c*0.8, c*1.2}}, {"d", 1, {0.5, 1.5}}};
}

std::shared_ptr<Fit> ExcludedVolumeFitter::fit() {
    fit_type = mini::type::DLIB_GLOBAL;
    setting::general::verbose = false;
    std::function<double(std::vector<double>)> f = std::bind(&ExcludedVolumeFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, f, guess);
    auto res = mini->minimize();

    update_excluded_volume(res.get_parameter("d").value);
    h.apply_water_scaling_factor(res.get_parameter("c").value);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    std::shared_ptr<Fit> ab_fit = fitter.fit();

    // update fitter object
    fitted = std::make_shared<Fit>(res, res.fval, data.size()-1); // start with the fit performed here
    fitted->add_fit(ab_fit);                                      // add the a,b inner fit
    fitted->add_plots(*this);                                     // make the result plottable
    fitted->evaluated_points = mini->get_evaluated_points();      // add the evaluated points

    setting::general::verbose = true;
    return fitted;
}

double ExcludedVolumeFitter::fit_only() {
    fit_type = mini::type::DLIB_GLOBAL;
    setting::general::verbose = false;
    std::function<double(std::vector<double>)> f = std::bind(&ExcludedVolumeFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, f, guess);
    auto res = mini->minimize();

    update_excluded_volume(res.get_parameter("d").value);
    h.apply_water_scaling_factor(res.get_parameter("c").value);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    return fitter.fit_only();
}

Fit::Plots ExcludedVolumeFitter::plot() {
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
    protein.update_effective_charge(d);
    h = protein.get_histogram();
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