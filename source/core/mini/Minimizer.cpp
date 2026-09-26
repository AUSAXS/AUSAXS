// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/Minimizer.h>

#include <mini/detail/Evaluation.h>
#include <mini/detail/Parameter.h>
#include <settings/GeneralSettings.h>
#include <utility/Exceptions.h>

#include <functional>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::mini;

Minimizer::Minimizer() = default;

Minimizer::~Minimizer() = default;

Minimizer::Minimizer(residual_function f) : objective(std::move(f)) {}

Result Minimizer::minimize() {
    if (!is_parameter_set()) {throw except::bad_order("Minimizer::minimize: No parameters were supplied.");}
    if (!is_function_set()) {throw except::bad_order("Minimizer::minimize: No function was set.");}

    clear_evaluated_points();
    return minimize_override();
}

void Minimizer::set_function(residual_function f) {
    objective = std::move(f);
}

std::vector<double> Minimizer::residuals(const std::vector<double>& params) {
    auto r = objective(params);
    fevals++;
    if (record) {evaluations.evals.emplace_back(params, chi2(r));}
    return r;
}

double Minimizer::function(const std::vector<double>& params) {
    return chi2(residuals(params));
}

double Minimizer::chi2(const std::vector<double>& r) {
    return std::transform_reduce(r.begin(), r.end(), 0.0, std::plus{}, [] (double v) {return v*v;});
}

Minimizer::residual_function Minimizer::get_recording_function() {
    return [this] (const std::vector<double>& params) {return residuals(params);};
}

bool Minimizer::empty() const noexcept {
    return !is_function_set() && !is_parameter_set();
}

void Minimizer::clear_parameters() noexcept {
    parameters.clear();
}

void Minimizer::record_evaluations(bool setting) {
    record = setting;
}

void Minimizer::add_parameter(const Parameter& param) {
    if (!param.has_bounds() && !param.has_guess()) {
        throw except::invalid_argument("Minimizer::add_parameter: Parameter \"" + param.name + "\"must either have a limit or a guess value.");
    }
    parameters.push_back(param);
}

void Minimizer::clear_evaluated_points() noexcept {
    evaluations.evals.clear();
}

bool Minimizer::is_function_set() const noexcept {
    return bool(objective); // functions are explicitly convertable to a bool which is true if a function has been set
}

bool Minimizer::is_parameter_set() const noexcept {
    return !parameters.empty();
}

mini::Landscape Minimizer::get_evaluated_points() const {
    if (evaluations.evals.empty()) {throw except::bad_order("Minimizer::get_evaluated_points: Cannot get evaluated points before a minimization call has been made.");}
    return evaluations;
}

mini::Landscape Minimizer::landscape(int bins) {
    if (parameters.empty()) {throw except::bad_order("Minimizer::landscape: No parameters were supplied.");}

    mini::Landscape l;
    const auto& bx_opt = parameters[0].bounds;
    if (!bx_opt.has_value()) {throw except::bad_order("Minimizer::landscape: Every scanned parameter must be bounded.");}
    auto bx = *bx_opt;
    for (int i = 0; i < bins; i++) {
        double vx = bx.min + i*bx.span()/(bins-1);
        double fval;
        if (parameters.size() == 2) {
            const auto& by_opt = parameters[1].bounds;
            if (!by_opt.has_value()) {throw except::bad_order("Minimizer::landscape: Every scanned parameter must be bounded.");}
            auto by = *by_opt;
            for (int j = 0; j < bins; j++) {
                double vy = by.min + j*by.span()/(bins-1);
                fval = function({vx, vy});
                l.evals.emplace_back(Evaluation{{vx, vy}, fval});
            }
        } else {
            fval = function({vx});
            l.evals.emplace_back(Evaluation{{vx}, fval});
        }

        // sanity check
        if (std::isnan(fval) || std::isinf(fval)) {
            if (settings::general::verbose) {std::cout << "Warning in Minimizer::landscape: Function value is nan or inf and will be skipped." << std::endl;}
            l.evals.pop_back();
            continue;
        }
    }
    return l;
}

void Minimizer::set_max_evals(int max_evals) {
    this->max_evals = max_evals;
}