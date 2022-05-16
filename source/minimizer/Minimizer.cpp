#include <minimizer/Minimizer.h>
#include <utility/Exceptions.h>

#include <functional>

using namespace mini;

Minimizer::Result::Result(const std::vector<double>& values, const std::vector<double>& errors, const std::vector<std::string>& names, double function_val) : fval(function_val) {
    if (values.size() != errors.size() || values.size() != names.size()) {throw except::size_error("Error in Minimizer::Result::Result: All provided vectors must be the same length.");}

    for (unsigned int i = 0; i < values.size(); i++) {
        this->params[names[i]] = values[i];
        this->errs[names[i]] = errors[i];
    }
}

Minimizer::Result::Result(const std::map<std::string, double>& params, const std::map<std::string, double>& errors, double function_val) : params(params), errs(errors), fval(function_val) {}

std::map<std::string, double> Minimizer::Result::parameters() const {return params;}

std::map<std::string, double> Minimizer::Result::errors() const {return errs;}

void Minimizer::set_function(double(&f)(double*)) {
    function = std::bind(f, std::placeholders::_1);
}

Minimizer::Minimizer(double(&f)(double*)) {
    set_function(f);
}