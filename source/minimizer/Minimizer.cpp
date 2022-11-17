#include <mini/Minimizer.h>
#include <utility/Exceptions.h>

#include <functional>

using namespace mini;

Minimizer::Minimizer(double(&f)(const double*)) {
    set_function(f);
}

Minimizer::Minimizer(std::function<double(const double*)> f) {
    set_function(f);
}

Result Minimizer::minimize() {
    if (!is_parameter_set()) {throw except::bad_order("Minimizer::minimize: No parameters were supplied.");}
    if (!is_function_set()) {throw except::bad_order("Minimizer::minimize: No function was set.");}

    clear_evaluated_points();
    return minimize_override();
}

void Minimizer::set_function(double(&f)(const double*)) {
    raw = std::bind(f, std::placeholders::_1);
    set_function(raw);
}

void Minimizer::set_function(std::function<double(const double*)> f) {
    raw = [f, this] (const double* par) {
        fevals++;
        return f(par);
    };
    wrapper = [this] (const double* par) {
        double fval = raw(par);
        std::vector<double> pars(par, par+parameters.size());
        evaluations.push_back(Evaluation(pars, fval));
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
    evaluations.clear();
}

bool Minimizer::is_function_set() const noexcept {
    return bool(function); // functions are explicitly convertable to a bool which is true if a function has been set
}

bool Minimizer::is_parameter_set() const noexcept {
    return !parameters.empty();
}