/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/SmartFitter.h>
#include <fitter/FitResult.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <math/CubicSpline.h>
#include <mini/All.h>
#include <settings/FitSettings.h>
#include <constants/ConstantsFitParameters.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::fitter;

unsigned int get_number_of_enabled_pars() {
    return 
        settings::fit::fit_hydration + 
        settings::fit::fit_excluded_volume + 
        settings::fit::fit_solvent_density + 
        settings::fit::fit_atomic_debye_waller + 
        settings::fit::fit_exv_debye_waller;
}

SmartFitter::~SmartFitter() = default;
SmartFitter::SmartFitter(SmartFitter&&) noexcept = default;
SmartFitter& SmartFitter::operator=(SmartFitter&&) noexcept = default;

SmartFitter::SmartFitter(const SimpleDataset& data) : data(data) {}

SmartFitter::SmartFitter(const SimpleDataset& saxs, std::unique_ptr<hist::DistanceHistogram> h) : SmartFitter(saxs) {
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
    if ((settings::fit::fit_atomic_debye_waller || settings::fit::fit_exv_debye_waller) && dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h) == nullptr) {
        throw except::invalid_argument("SmartFitter::validate_model: Histogram must support Debye-Waller operations.");
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
    if (settings::fit::fit_hydration) {
        guess.push_back(mini::Parameter{constants::fit::to_string(constants::fit::Parameters::SCALING_WATER), 1, cast_h(model.get())->get_water_scaling_factor_limits()});
    }

    if (settings::fit::fit_excluded_volume) {
        guess.push_back(mini::Parameter{constants::fit::to_string(constants::fit::Parameters::SCALING_EXV), 1, cast_exv(model.get())->get_excluded_volume_scaling_factor_limits()});
    }

    if (settings::fit::fit_solvent_density) {
        guess.push_back(mini::Parameter{constants::fit::to_string(constants::fit::Parameters::SCALING_RHO), 1, cast_exv(model.get())->get_solvent_density_scaling_factor_limits()});
    }

    if (settings::fit::fit_atomic_debye_waller) {
        guess.push_back(mini::Parameter{constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_ATOMIC), 0, cast_exv(model.get())->get_debye_waller_factor_limits()});
    }

    if (settings::fit::fit_exv_debye_waller) {
        guess.push_back(mini::Parameter{constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_EXV), 0, cast_exv(model.get())->get_debye_waller_factor_limits()});
    }
    return guess;
}

fitter::detail::LinearLeastSquares SmartFitter::prepare_linear_fitter(const std::vector<double>& params) {
    assert(
        params.size() == get_number_of_enabled_pars()
        && "SmartFitter::get_model_curve: Invalid number of parameters."
    );

    int index = 0;
    if (settings::fit::fit_hydration)           {cast_h(model.get())->apply_water_scaling_factor(params[index++]);}
    if (settings::fit::fit_excluded_volume)     {cast_exv(model.get())->apply_excluded_volume_scaling_factor(params[index++]);}
    if (settings::fit::fit_solvent_density)     {cast_exv(model.get())->apply_solvent_density_scaling_factor(params[index++]);}
    if (settings::fit::fit_atomic_debye_waller) {cast_exv(model.get())->apply_atomic_debye_waller_factor(params[index++]);}
    if (settings::fit::fit_exv_debye_waller)    {cast_exv(model.get())->apply_exv_debye_waller_factor(params[index++]);}
    return detail::LinearLeastSquares(splice(model->debye_transform().get_counts()), data.y(), data.yerr());
}

std::unique_ptr<FitResult> SmartFitter::fit() {
    validate_model(model.get());
    if (guess.empty()) {guess = get_default_guess();}

    auto f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(algorithm, std::move(f), guess);
    auto res = mini->minimize();

    auto linear_fitter = prepare_linear_fitter(res.get_parameter_values());
    auto linear_fit = linear_fitter.fit();

    assert(std::abs(linear_fit->fval - res.fval) < 1e-6 && "SmartFitter::fit: Linear fit and minimizer results do not match.");
    auto fit_result = std::make_unique<FitResult>(res, res.fval, dof()+2);     // start with the fit performed here
    fit_result->add_fit(linear_fit.get(), true);                               // add the a,b inner fit
    fit_result->set_data_curves(
        data.x(),                                                              // q
        data.y(),                                                              // I
        data.yerr(),                                                           // I_err
        linear_fitter.get_model_curve(linear_fit->get_parameter_values()),     // I_fit
        linear_fitter.get_residuals(linear_fit->get_parameter_values())        // residuals
    );
    fit_result->evaluated_points = mini->get_evaluated_points();               // add the evaluated points
    return fit_result;
}

std::vector<double> SmartFitter::fit_params_only() {
    validate_model(model.get());
    if (guess.empty()) {guess = get_default_guess();}

    std::function<double(std::vector<double>)> f = std::bind(&SmartFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(algorithm, std::move(f), guess);
    return mini->minimize().get_parameter_values();
}

std::vector<double> SmartFitter::get_residuals(const std::vector<double>& params) {
    auto fitter = prepare_linear_fitter(params);
    return fitter.get_residuals();
}

std::vector<double> SmartFitter::get_model_curve(const std::vector<double>& params) {
    auto fitter = prepare_linear_fitter(params);
    return fitter.get_model_curve();
}

observer_ptr<hist::DistanceHistogram> SmartFitter::get_model() {
    return model.get();
}

SimpleDataset SmartFitter::get_data() const {
    return data;
}

unsigned int SmartFitter::size() const {
    return data.size();
}

unsigned int SmartFitter::dof() const {
    return data.size() - 2 - get_number_of_enabled_pars();
}

void SmartFitter::set_guess(std::vector<mini::Parameter>&& guess) {
    if (unsigned int N = get_number_of_enabled_pars(); guess.size() != N) {
        throw except::invalid_argument("SmartFitter::set_guess: Invalid number of parameters. Got " + std::to_string(guess.size()) + ", expected " + std::to_string(N) + ".");
    }

    // validate and reorder the parameters
    std::vector<int> order;
    for (unsigned int i = 0; i < guess.size(); ++i) {
        if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::SCALING_WATER)) {
            if (!settings::fit::fit_hydration) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set hydration scaling factor when hydration is disabled.");}
            order.push_back(0);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::SCALING_EXV)) {
            if (!settings::fit::fit_excluded_volume) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set excluded volume scaling factor when excluded volume is disabled.");}
            order.push_back(1);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::SCALING_RHO)) {
            if (!settings::fit::fit_solvent_density) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set solvent density scaling factor when solvent density is disabled.");}
            order.push_back(2);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_ATOMIC)) {
            if (!settings::fit::fit_atomic_debye_waller) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set atomic Debye-Waller factor when atomic Debye-Waller is disabled.");}
            order.push_back(3);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_EXV)) {
            if (!settings::fit::fit_exv_debye_waller) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set excluded volume Debye-Waller factor when excluded volume Debye-Waller is disabled.");}
            order.push_back(4);
        } else {
            throw except::invalid_argument("SmartFitter::set_guess: Unknown parameter name: \"" + guess[i].name + "\"");
        }
    }
    this->guess.clear();
    std::for_each(order.begin(), order.end(), [&] (int i) {this->guess.push_back(guess[i]);});
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