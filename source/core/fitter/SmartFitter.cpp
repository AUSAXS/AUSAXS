/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/SmartFitter.h>
#include <fitter/LinearFitter.h>
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
    if (settings::fit::fit_solvent_density) {guess.push_back(mini::Parameter{"e", 1, cast_exv(model.get())->get_solvent_density_scaling_factor_limits()});}
    return guess;
}

std::vector<double> SmartFitter::extract_opt_pars(observer_ptr<FitResult> smart) {
    std::vector<double> popt;
    if (settings::fit::fit_hydration) {popt.push_back(smart->get_parameter("c"));}
    if (settings::fit::fit_excluded_volume) {popt.push_back(smart->get_parameter("d"));}
    if (settings::fit::fit_solvent_density) {popt.push_back(smart->get_parameter("e"));}
    return popt;
};

std::unique_ptr<FitResult> SmartFitter::fit() {
    validate_model(model.get());
    if (guess.empty()) {guess = get_default_guess();}

    auto f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(algorithm, std::move(f), guess);
    auto res = mini->minimize();

    auto ym = get_model_curve(res.get_parameter_values());
    LinearFitter fitter(ym, data.y(), data.yerr());
    std::shared_ptr<FitResult> ab_fit = fitter.fit();

    auto fit_result = std::make_unique<FitResult>(res, res.fval, data.size()-1); // start with the fit performed here
    fit_result->add_fit(ab_fit.get(), true);                                     // add the a,b inner fit
    fit_result->set_data_curves(
        data.x(),                                                                // q
        data.y(),                                                                // I
        data.yerr(),                                                             // I_err
        std::move(ym),                                                           // I_fit
        get_residuals(extract_opt_pars(fit_result.get()))                        // residuals
    );
    fit_result->evaluated_points = mini->get_evaluated_points();                 // add the evaluated points
    return fit_result;
}

double SmartFitter::fit_chi2_only() {
    validate_model(model.get());
    if (guess.empty()) {guess = get_default_guess();}

    std::function<double(std::vector<double>)> f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(algorithm, std::move(f), guess);
    return mini->minimize().fval;
}

std::vector<double> SmartFitter::get_residuals(const std::vector<double>& params) {
    auto Im = get_model_curve(params);
    LinearFitter fitter(Im, data.y(), data.yerr());
    return fitter.get_residuals();
}

std::vector<double> SmartFitter::get_model_curve(const std::vector<double>& params) {
    assert(
        static_cast<int>(params.size()) == settings::fit::fit_hydration + settings::fit::fit_excluded_volume + settings::fit::fit_solvent_density 
        && "SmartFitter::get_model_curve: Invalid number of parameters."
    );

    int index = 0;
    if (settings::fit::fit_hydration) {cast_h(model.get())->apply_water_scaling_factor(params[index++]);}
    if (settings::fit::fit_excluded_volume) {cast_exv(model.get())->apply_excluded_volume_scaling_factor(params[index++]);}
    if (settings::fit::fit_solvent_density) {cast_exv(model.get())->apply_solvent_density_scaling_factor(params[index]);}

    return splice(model->debye_transform().get_counts());
}

SimpleDataset SmartFitter::get_data() const {
    return data;
}

unsigned int SmartFitter::size() const {
    return data.size();
}

unsigned int SmartFitter::dof() const {
    return data.size() - 2 - settings::fit::fit_excluded_volume - settings::fit::fit_hydration - settings::fit::fit_solvent_density;
}

void SmartFitter::set_guess(std::vector<mini::Parameter>&& guess) {
    if (unsigned int N = settings::fit::fit_hydration + settings::fit::fit_excluded_volume + settings::fit::fit_solvent_density; guess.size() != N) {
        throw except::invalid_argument("SmartFitter::set_guess: Invalid number of parameters. Got " + std::to_string(guess.size()) + ", expected " + std::to_string(N) + ".");
    }
    for (auto g : guess) {
        if (g.name == "c" && !settings::fit::fit_hydration) {
            throw except::invalid_argument("SmartFitter::set_guess: Cannot set hydration scaling factor when hydration is disabled.");
        }
        if (g.name == "d" && !settings::fit::fit_excluded_volume) {
            throw except::invalid_argument("SmartFitter::set_guess: Cannot set excluded volume scaling factor when excluded volume is disabled.");
        }
        if (g.name == "e" && !settings::fit::fit_solvent_density) {
            throw except::invalid_argument("SmartFitter::set_guess: Cannot set solvent density scaling factor when solvent density is disabled.");
        }
    }
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