/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/refactored/SmartFitter.h>
#include <fitter/refactored/LinearFitter.h>
#include <fitter/FitResult.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <math/CubicSpline.h>
#include <mini/All.h>
#include <settings/FitSettings.h>

#include <cassert>

using namespace fitter;

SmartFitter::~SmartFitter() = default;

SmartFitter::SmartFitter(const SimpleDataset& data) : data(data) {}

SmartFitter::SmartFitter(const SimpleDataset& saxs, std::unique_ptr<hist::ICompositeDistanceHistogram> h) : SmartFitter(saxs) {
    set_model(std::move(h));
}

void validate_model(observer_ptr<hist::DistanceHistogram> h) {
    if (h == nullptr) {throw except::invalid_argument("SmartFitter::validate_model: Cannot fit without a model.");}
    if (settings::fit::fit_hydration && dynamic_cast<hist::ICompositeDistanceHistogram*>(h) == nullptr) {
        throw except::invalid_argument("SmartFitter::validate_model: Histogram must support hydration shell operations.");
    }
    if ((settings::fit::fit_excluded_volume || settings::fit::fit_solvent_density) && dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h) == nullptr) {
        throw except::invalid_argument("SmartFitter::validate_model: Histogram must support excluded volume operations.");
    }
}

observer_ptr<hist::ICompositeDistanceHistogramExv> cast_exv(observer_ptr<hist::DistanceHistogram> hist) {
    return static_cast<hist::ICompositeDistanceHistogramExv*>(hist);
}

observer_ptr<hist::ICompositeDistanceHistogram> cast_h(observer_ptr<hist::DistanceHistogram> hist) {
    return static_cast<hist::ICompositeDistanceHistogram*>(hist);
}

std::vector<mini::Parameter> SmartFitter::get_default_guess() const {
    std::vector<mini::Parameter> guess;
    if (settings::fit::fit_hydration) {guess.push_back(mini::Parameter{"c", 1, cast_h(model.get())->get_water_scaling_factor_limits()});}
    if (settings::fit::fit_excluded_volume) {guess.push_back(mini::Parameter{"d", 1, cast_exv(model.get())->get_excluded_volume_scaling_factor_limits()});}
    if (settings::fit::fit_solvent_density) {guess.push_back(mini::Parameter{"e", 1, {0.95, 1.05}});}
    return guess;
}

std::vector<double> SmartFitter::extract_opt_pars(observer_ptr<FitResult> smart) {
    std::vector<double> popt = {smart->get_parameter(0), smart->get_parameter(1)};
    if (settings::fit::fit_hydration) {popt.push_back(smart->get_parameter("c"));}
    if (settings::fit::fit_excluded_volume) {popt.push_back(smart->get_parameter("d"));}
    if (settings::fit::fit_solvent_density) {popt.push_back(smart->get_parameter("e"));}
    return popt;
};

std::unique_ptr<FitResult> SmartFitter::fit() {
    validate_model(model.get());

    // we interpolate the data to avoid having to interpolate the model for every chi2 evaluation
    // for the final chi2 evaluation we use the original data
    data_spliced = data.interpolate(model->get_q_axis());
    data_spliced.yerr() = data.yerr(); // keep the errors since we have no guarantee of smoothness

    auto f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(mini::type::DEFAULT, std::move(f), guess);
    auto res = mini->minimize();

    if (settings::fit::fit_hydration) {cast_h(model.get())->apply_water_scaling_factor(res.get_parameter("c"));}
    if (settings::fit::fit_excluded_volume) {cast_exv(model.get())->apply_excluded_volume_scaling_factor(res.get_parameter("d"));}
    if (settings::fit::fit_solvent_density) {cast_exv(model.get())->apply_solvent_density_scaling_factor(res.get_parameter("e"));}

    auto y_model_i = splice(model->debye_transform().get_counts());
    LinearFitter fitter(y_model_i, data.y(), data.yerr());
    std::shared_ptr<FitResult> ab_fit = fitter.fit();

    auto fit_result = std::make_unique<FitResult>(res, res.fval, data.size()-1); // start with the fit performed here
    fit_result->add_fit(ab_fit.get());                                           // add the a,b inner fit
    fit_result->curves = {{data.x(), data.y(), y_model_i, get_residuals(extract_opt_pars(fit_result.get()))}};
    fit_result->evaluated_points = mini->get_evaluated_points();                 // add the evaluated points
    return fit_result;
}

double SmartFitter::fit_chi2_only() {
    validate_model(model.get());

    std::function<double(std::vector<double>)> f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(mini::type::DEFAULT, std::move(f), guess);
    auto res = mini->minimize();

    if (settings::fit::fit_hydration) {cast_h(model.get())->apply_water_scaling_factor(res.get_parameter("c"));}
    if (settings::fit::fit_excluded_volume) {cast_exv(model.get())->apply_excluded_volume_scaling_factor(res.get_parameter("d"));}
    if (settings::fit::fit_solvent_density) {cast_exv(model.get())->apply_solvent_density_scaling_factor(res.get_parameter("e"));}

    LinearFitter fitter(splice(model->debye_transform().get_counts()), data.y(), data.yerr());
    return fitter.fit_chi2_only();
}

std::vector<double> SmartFitter::get_residuals(const std::vector<double>& params) {
    assert(params.size() == settings::fit::fit_hydration + settings::fit::fit_excluded_volume + settings::fit::fit_solvent_density && "SmartFitter::chi2: Invalid number of parameters.");

    int index = 0;
    if (settings::fit::fit_hydration) {cast_h(model.get())->apply_water_scaling_factor(params[index++]);}
    if (settings::fit::fit_excluded_volume) {cast_exv(model.get())->apply_excluded_volume_scaling_factor(params[index++]);}
    if (settings::fit::fit_solvent_density) {cast_exv(model.get())->apply_solvent_density_scaling_factor(params[index]);}

    LinearFitter fitter(model->debye_transform().get_counts(), data_spliced.y(), data_spliced.yerr());
    return fitter.get_residuals(params);
}

SimpleDataset SmartFitter::get_data() const {
    return data;
}

unsigned int SmartFitter::dof() const {
    return data.size() - 2 - settings::fit::fit_excluded_volume - settings::fit::fit_hydration - settings::fit::fit_solvent_density;
}

void SmartFitter::set_guess(std::vector<mini::Parameter>&& guess) {
    this->guess = std::move(guess);
}

void SmartFitter::set_model(std::unique_ptr<hist::DistanceHistogram> h) {
    model = std::move(h);
}

std::vector<double> SmartFitter::splice(const std::vector<double>& ym) const {
    std::vector<double> Im(data.size()); // spliced model values
    math::CubicSpline s(model->get_q_axis(), ym);
    for (unsigned int i = 0; i < data.size(); ++i) {
        Im[i] = s.spline(data.x(i));
    }
    return Im;
}