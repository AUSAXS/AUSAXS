/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/refactored/SmartFitter.h>
#include <fitter/refactored/SimpleLeastSquares.h>
#include <fitter/FitResult.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <mini/All.h>

#include <cassert>

using namespace fitter;

SmartFitter::~SmartFitter() = default;

SmartFitter::SmartFitter(const SimpleDataset& data) : data(data) {}

SmartFitter::SmartFitter(const SimpleDataset& saxs, std::unique_ptr<hist::ICompositeDistanceHistogram> h) : SmartFitter(saxs) {
    set_model(std::move(h));
}

bool fit_exv, fit_hydration, fit_solvent_density;
void validate_model(observer_ptr<hist::ICompositeDistanceHistogram> h) {
    if (h == nullptr) {throw except::invalid_argument("SmartFitter::validate_model: Cannot fit without a model.");}
    if (fit_exv && dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h) == nullptr) {
        throw except::invalid_argument("SmartFitter::validate_model: Histogram must support excluded volume operations.");
    }
}

observer_ptr<hist::ICompositeDistanceHistogramExv> cast_exv(observer_ptr<hist::ICompositeDistanceHistogram> hist) {
    return static_cast<hist::ICompositeDistanceHistogramExv*>(hist);
}

observer_ptr<hist::ICompositeDistanceHistogram> cast_h(observer_ptr<hist::DistanceHistogram> hist) {
    return static_cast<hist::ICompositeDistanceHistogram*>(hist);
}

std::vector<mini::Parameter> SmartFitter::get_default_guess() const {
    std::vector<mini::Parameter> guess;
    if (fit_hydration) {guess.push_back(mini::Parameter{"c", 1, model->get_water_scaling_factor_limits()});}
    if (fit_exv) {guess.push_back(mini::Parameter{"d", 1, cast_exv(model.get())->get_excluded_volume_scaling_factor_limits()});}
    if (fit_solvent_density) {guess.push_back(mini::Parameter{"e", 1, {0.95, 1.05}});}
    return guess;
}

SimpleDataset SmartFitter::get_optimized_model() const {}

std::unique_ptr<FitResult> SmartFitter::fit() {
    validate_model(model.get());

    auto f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(mini::type::DEFAULT, std::move(f), guess);
    auto res = mini->minimize();

    if (fit_hydration) {cast_h(model.get())->apply_water_scaling_factor(res.get_parameter("c"));}
    if (fit_exv) {cast_exv(model.get())->apply_excluded_volume_scaling_factor(res.get_parameter("d"));}
    if (fit_solvent_density) {(model.get())->apply_solvent_density_scaling_factor(res.get_parameter("e"));}

    std::vector<double> ym = model->debye_transform().get_counts();
    SimpleLeastSquares fitter(splice(ym), data.y(), data.yerr());
    std::shared_ptr<FitResult> ab_fit = fitter.fit();

    // update fitter object
    auto fit_result = std::make_unique<FitResult>(res, res.fval, data.size()-1); // start with the fit performed here
    fit_result->add_fit(ab_fit.get());                                           // add the a,b inner fit
    fit_result->add_plots(this);                                                 // make the result plottable
    fit_result->evaluated_points = mini->get_evaluated_points();                 // add the evaluated points
    return fit_result;
}

double SmartFitter::fit_chi2_only() {
    validate_model(model.get());

    std::function<double(std::vector<double>)> f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(mini::type::DEFAULT, std::move(f), guess);
    auto res = mini->minimize();

    if (fit_hydration) {cast_h(model.get())->apply_water_scaling_factor(res.get_parameter("c"));}
    if (fit_exv) {cast_exv(model.get())->apply_excluded_volume_scaling_factor(res.get_parameter("d"));}
    if (fit_solvent_density) {(model.get())->apply_solvent_density_scaling_factor(res.get_parameter("e"));}

    std::vector<double> ym = model->debye_transform().get_counts();
    SimpleLeastSquares fitter(splice(ym), data.y(), data.yerr());
    return fitter.fit_chi2_only();
}

double SmartFitter::chi2(const std::vector<double>& params) {
    assert(params.size() == fit_hydration + fit_exv + fit_solvent_density && "SmartFitter::chi2: Invalid number of parameters.");

    int index = 0;
    if (fit_hydration) {model->apply_water_scaling_factor(params[index++]);}
    if (fit_exv) {cast_exv(model.get())->apply_excluded_volume_scaling_factor(params[index++]);}
    if (fit_solvent_density) {model->apply_solvent_density_scaling_factor(params[index]);}

    update_excluded_volume(params[1]);
    return HydrationFitter::chi2(params);
}

SimpleDataset SmartFitter::get_data() const {
    return data;
}

unsigned int SmartFitter::dof() const {
    return data.size() - 2 - fit_exv - fit_hydration - fit_solvent_density;
}

void SmartFitter::set_guess(std::vector<mini::Parameter>&& guess) {
    this->guess = std::move(guess);
}

void SmartFitter::set_model(std::unique_ptr<hist::DistanceHistogram> h) {
    HydrationFitter hfit(data, std::move(h));
    auto hres = hfit.fit();
    double c = hres->get_parameter("c").value;
    HydrationFitter::operator=(std::move(hfit));
    initialize_guess();
    validate_histogram();
    guess.guess = c;
}