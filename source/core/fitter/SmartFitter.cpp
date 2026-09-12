// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <fitter/SmartFitter.h>

#include <constants/ConstantsFitParameters.h>
#include <dataset/SimpleDataset.h>
#include <fitter/FitResult.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <math/CubicSpline.h>
#include <mini/MinimizerFactory.h>
#include <settings/FitSettings.h>
#include <utility/Console.h>

#include <algorithm>
#include <cassert>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::fitter;

namespace {
    void warn_if_parameter_on_bound(const std::vector<mini::Parameter>& guess, const mini::Result& res) {
        constexpr double rel_tol = 1e-3;    // fraction of the allowed range counted as "on the bound"
        for (int i = 0; i < static_cast<int>(guess.size()); ++i) {
            const auto& param = guess[i];
            if (!param.bounds.has_value()) {continue;}
            const auto& bounds = param.bounds.value();
            double span = bounds.span();
            if (!(0 < span)) {continue;}
            double value = res.get_parameter(i).value;
            double tol = rel_tol*span;
            bool on_lower = value <= bounds.min + tol;
            bool on_upper = bounds.max - tol <= value;
            if (!on_lower && !on_upper) {continue;}
            console::print_warning(
                "Warning: the fitted parameter \"" + param.name + "\" converged to its " + (on_lower ? "lower" : "upper") + 
                " bound (" + std::to_string(value) + " in [" + std::to_string(bounds.min) + ", " + std::to_string(bounds.max) + "])."
            );
        }
    }
    bool solvent_density_warned = false;
}

SmartFitter::~SmartFitter() {solvent_density_warned = false;}
SmartFitter::SmartFitter(SmartFitter&&) noexcept = default;
SmartFitter& SmartFitter::operator=(SmartFitter&&) noexcept = default;

SmartFitter::EnabledFitParameters SmartFitter::EnabledFitParameters::initialize_from_settings() {
    return {
        .hydration = settings::fit::fit_hydration,
        .excluded_volume = settings::fit::fit_excluded_volume,
        .solvent_density = settings::fit::fit_solvent_density,
        .atomic_debye_waller = settings::fit::fit_atomic_debye_waller,
        .exv_debye_waller = settings::fit::fit_exv_debye_waller
    };
}

int SmartFitter::EnabledFitParameters::get_enabled_pars_count() const {return static_cast<int>(hydration)+static_cast<int>(excluded_volume)+static_cast<int>(solvent_density)+static_cast<int>(atomic_debye_waller)+static_cast<int>(exv_debye_waller);}

void SmartFitter::EnabledFitParameters::apply_pars(const std::vector<double>& params, observer_ptr<hist::DistanceHistogram> model) const {
    assert(
        static_cast<int>(params.size()) == get_enabled_pars_count()
        && "SmartFitter::EnabledFitParameters::apply_pars: Invalid number of parameters."
    );

    int index = 0;
    if (hydration)           {static_cast<hist::ICompositeDistanceHistogram*>(model)->apply_water_scaling_factor(params[index++]);}
    if (excluded_volume)     {static_cast<hist::ICompositeDistanceHistogramExv*>(model)->apply_excluded_volume_scaling_factor(params[index++]);}
    if (solvent_density)     {static_cast<hist::ICompositeDistanceHistogramExv*>(model)->apply_solvent_density_scaling_factor(params[index++]);}
    if (atomic_debye_waller) {static_cast<hist::ICompositeDistanceHistogramExv*>(model)->apply_atomic_debye_waller_factor(params[index++]);}
    if (exv_debye_waller)    {static_cast<hist::ICompositeDistanceHistogramExv*>(model)->apply_exv_debye_waller_factor(params[index++]);}
}

void SmartFitter::EnabledFitParameters::validate_model(observer_ptr<hist::DistanceHistogram> h) {
    if (h == nullptr) {throw except::invalid_argument("SmartFitter::EnabledFitParameters::validate_model: Cannot fit without a model.");}
    if (hydration && dynamic_cast<hist::ICompositeDistanceHistogram*>(h) == nullptr) {
        console::print_warning("SmartFitter::EnabledFitParameters::validate_model: Hydration shell fitting is enabled, but the model does not support hydration shell operations. Disabling hydration shell fitting.");
        hydration = false;
    }
    if ((excluded_volume || solvent_density) && dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h) == nullptr) {
        console::print_warning("SmartFitter::EnabledFitParameters::validate_model: Excluded volume fitting is enabled, but the model does not support excluded volume operations. Disabling excluded volume fitting.");
        excluded_volume = false;
        solvent_density = false;
    }
    if ((atomic_debye_waller || exv_debye_waller) && dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h) == nullptr) {
        console::print_warning("SmartFitter::EnabledFitParameters::validate_model: Debye-Waller fitting is enabled, but the model does not support Debye-Waller operations. Disabling Debye-Waller fitting.");
        atomic_debye_waller = false;
        exv_debye_waller = false;
    }
    if (solvent_density) {
        if (!solvent_density_warned) {
            console::print_warning(
                "Warning: fitting the solvent density scaling factor "
                "(\"" + constants::fit::to_string(constants::fit::Parameters::SCALING_RHO) + "\") is strongly discouraged."
            );
            solvent_density_warned = true;
        }
    }
}

SmartFitter::SmartFitter(SimpleDataset data) : data(std::move(data)) {
    enabled_fit_parameters = EnabledFitParameters::initialize_from_settings();
    solvent_density_warned = false;
}

SmartFitter::SmartFitter(const SimpleDataset& data, std::unique_ptr<hist::DistanceHistogram> h) : SmartFitter(data) {
    enabled_fit_parameters = EnabledFitParameters::initialize_from_settings();
    set_model(std::move(h));
    solvent_density_warned = false;
}

namespace {
    observer_ptr<hist::ICompositeDistanceHistogramExv> cast_exv(observer_ptr<hist::DistanceHistogram> hist) {
        return static_cast<hist::ICompositeDistanceHistogramExv*>(hist);
    }

    observer_ptr<hist::ICompositeDistanceHistogram> cast_h(observer_ptr<hist::DistanceHistogram> hist) {
        return static_cast<hist::ICompositeDistanceHistogram*>(hist);
    }
}

std::vector<mini::Parameter> SmartFitter::get_default_guess() const {
    std::vector<mini::Parameter> guess;
    if (enabled_fit_parameters.hydration) {
        guess.emplace_back(
            constants::fit::to_string(constants::fit::Parameters::SCALING_WATER), 
            1, 
            cast_h(model.get())->get_water_scaling_factor_limits()
        );
    }

    if (enabled_fit_parameters.excluded_volume) {
        guess.emplace_back(
            constants::fit::to_string(constants::fit::Parameters::SCALING_EXV), 
            1, 
            cast_exv(model.get())->get_excluded_volume_scaling_factor_limits()
        );
    }

    if (enabled_fit_parameters.solvent_density) {
        guess.emplace_back(
            constants::fit::to_string(constants::fit::Parameters::SCALING_RHO), 
            1, 
            cast_exv(model.get())->get_solvent_density_scaling_factor_limits()
        );
    }

    if (enabled_fit_parameters.atomic_debye_waller) {
        guess.emplace_back(
            constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_ATOMIC), 
            0, 
            cast_exv(model.get())->get_debye_waller_factor_limits()
        );
    }

    if (enabled_fit_parameters.exv_debye_waller) {
        guess.emplace_back(
            constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_EXV), 
            0, 
            cast_exv(model.get())->get_debye_waller_factor_limits()
        );
    }
    return guess;
}

fitter::detail::LinearLeastSquares SmartFitter::prepare_linear_fitter(const std::vector<double>& params) {
    assert(
        static_cast<int>(params.size()) == enabled_fit_parameters.get_enabled_pars_count()
        && "SmartFitter::get_model_curve: Invalid number of parameters."
    );
    enabled_fit_parameters.apply_pars(params, model.get());
    return {splice(model->debye_transform().get_counts()), data.y(), data.yerr()};
}

std::unique_ptr<FitResult> SmartFitter::fit() {
    enabled_fit_parameters = EnabledFitParameters::initialize_from_settings();
    enabled_fit_parameters.validate_model(model.get());
    if (guess.empty()) {guess = get_default_guess();}

    if (enabled_fit_parameters.get_enabled_pars_count() == 0) {
        auto linear_fitter = prepare_linear_fitter({});
        return linear_fitter.fit();    
    }

    auto f = [this] (const std::vector<double>& params) {return chi2(params);};
    auto mini = mini::create_minimizer(algorithm, std::move(f), guess);
    auto res = mini->minimize();
    warn_if_parameter_on_bound(guess, res);

    auto linear_fitter = prepare_linear_fitter(res.get_parameter_values());
    auto linear_fit = linear_fitter.fit();

    auto fit_result = std::make_unique<FitResult>(res, dof()+2);            // start with the fit performed here
    fit_result->add_fit(linear_fit.get(), true);                            // add the a,b inner fit
    fit_result->set_data_curves(
        data.x(),                                                           // q
        data.y(),                                                           // I
        data.yerr(),                                                        // I_err
        linear_fitter.get_model_curve(linear_fit->get_parameter_values()),  // I_fit
        linear_fitter.get_residuals(linear_fit->get_parameter_values())     // residuals
    );
    fit_result->evaluated_points = mini->get_evaluated_points();            // add the evaluated points
    return fit_result;
}

std::vector<double> SmartFitter::fit_params_only() {
    enabled_fit_parameters = EnabledFitParameters::initialize_from_settings();
    enabled_fit_parameters.validate_model(model.get());
    if (guess.empty()) {guess = get_default_guess();}

    std::function<double(std::vector<double>)> f = [this](auto && PH1) { return chi2(std::forward<decltype(PH1)>(PH1)); };
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

int SmartFitter::size() const {
    return data.size();
}

int SmartFitter::dof() const {
    return data.size() - 2 - enabled_fit_parameters.get_enabled_pars_count();
}

void SmartFitter::set_guess(const std::vector<mini::Parameter>& guess) {
    if (int N = enabled_fit_parameters.get_enabled_pars_count(); static_cast<int>(guess.size()) != N) {
        throw except::invalid_argument("SmartFitter::set_guess: Invalid number of parameters. Got " + std::to_string(guess.size()) + ", expected " + std::to_string(N) + ".");
    }

    // validate and reorder the parameters
    // note: 'order' pairs the canonical slot of each supplied parameter with its index in the input,
    //       so that sorting it yields the permutation taking the input into canonical order
    std::vector<std::pair<int, int>> order;
    for (int i = 0; i < static_cast<int>(guess.size()); ++i) {
        if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::SCALING_WATER)) {
            if (!enabled_fit_parameters.hydration) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set hydration scaling factor when hydration is disabled.");}
            order.emplace_back(0, i);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::SCALING_EXV)) {
            if (!enabled_fit_parameters.excluded_volume) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set excluded volume scaling factor when excluded volume is disabled.");}
            order.emplace_back(1, i);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::SCALING_RHO)) {
            if (!enabled_fit_parameters.solvent_density) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set solvent density scaling factor when solvent density is disabled.");}
            order.emplace_back(2, i);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_ATOMIC)) {
            if (!enabled_fit_parameters.atomic_debye_waller) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set atomic Debye-Waller factor when atomic Debye-Waller is disabled.");}
            order.emplace_back(3, i);
        } else if (guess[i].name == constants::fit::to_string(constants::fit::Parameters::DEBYE_WALLER_EXV)) {
            if (!enabled_fit_parameters.exv_debye_waller) {throw except::invalid_argument("SmartFitter::set_guess: Cannot set excluded volume Debye-Waller factor when excluded volume Debye-Waller is disabled.");}
            order.emplace_back(4, i);
        } else {
            throw except::invalid_argument("SmartFitter::set_guess: Unknown parameter name: \"" + guess[i].name + "\"");
        }
    }
    std::ranges::sort(order);
    assert(
        std::ranges::adjacent_find(order, [] (const auto& a, const auto& b) {return a.first == b.first;}) == order.end()
        && "SmartFitter::set_guess: The same parameter was supplied more than once."
    );

    this->guess.clear();
    this->guess.reserve(order.size());
    std::ranges::for_each(order, [&] (const auto& o) {this->guess.push_back(std::move(guess[o.second]));});
}

void SmartFitter::set_model(std::unique_ptr<hist::DistanceHistogram> h) {
    model = std::move(h);
}

std::vector<double> SmartFitter::splice(const std::vector<double>& ym) const {
    std::vector<double> Im(data.size()); // spliced model values
    math::CubicSpline s(hist::DistanceHistogram::get_q_axis(), ym);
    for (int i = 0; i < data.size(); ++i) {
        Im[i] = s.spline(data.x(i));
    }
    return Im;
}