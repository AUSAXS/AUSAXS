// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/Minimizer.h>

#include <mini/detail/Evaluation.h>
#include <mini/detail/Parameter.h>
#include <settings/GeneralSettings.h>
#include <utility/Exceptions.h>

#include <functional>

using namespace ausaxs;
using namespace ausaxs::mini;

namespace {
    // a plain function has no state to own, so it is captured as a pointer rather than by reference
    std::function<double(std::vector<double>)> as_function(double(&f)(std::vector<double>)) {
        return [fp = &f] (std::vector<double> p) {return fp(std::move(p));};
    }
}

Minimizer::Minimizer() = default;

Minimizer::~Minimizer() = default;

Minimizer::Minimizer(double(&f)(std::vector<double>)) {
    _set_function(as_function(f));
}

Minimizer::Minimizer(std::function<double(std::vector<double>)>&& f) {
    _set_function(std::move(f));
}

Result Minimizer::minimize() {
    if (!is_parameter_set()) {throw except::bad_order("Minimizer::minimize: No parameters were supplied.");}
    if (!is_function_set()) {throw except::bad_order("Minimizer::minimize: No function was set.");}

    clear_evaluated_points();
    return minimize_override();
}

void Minimizer::set_function(double(&f)(std::vector<double>)) {
    set_function(as_function(f));
}

void Minimizer::set_function(std::function<double(std::vector<double>)>&& f) {
    _set_function(std::move(f));
}

void Minimizer::_set_function(std::function<double(std::vector<double>)>&& f) {
    raw = std::move(f);
    wrapper = [this] (std::vector<double> p) {
        double fval = raw(p);
        evaluations.evals.emplace_back(std::move(p), fval);
        fevals++;
        return fval;
    };

    function = wrapper;
}

bool Minimizer::empty() const noexcept {
    return !is_function_set() && !is_parameter_set();
}

void Minimizer::clear_parameters() noexcept {
    parameters.clear();
}

void Minimizer::record_evaluations(bool setting) {
    function = setting ? wrapper : raw;
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
    return bool(function); // functions are explicitly convertable to a bool which is true if a function has been set
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