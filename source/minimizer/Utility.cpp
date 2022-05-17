#include <minimizer/Utility.h>

using namespace mini;

Result::Result(const FittedParameter& param, double fval) : parameters({param}), fval(fval) {}
Result::Result(const std::vector<FittedParameter>& params, double fval) : parameters(params), fval(fval) {}
Result::get_parameter(std::string name);


Evaluation::Evaluation(std::vector<double> vals, double fval) : vals(vals), fval(fval) {}


Parameter::Parameter(std::string name, double guess = 0, Limit bounds = {0, 0}) : name(name), guess(guess), bounds(bounds) {}
bool Parameter::has_bounds() const noexcept {return bounds.has_value();}
bool Parameter::has_guess() const noexcept {return guess.has_value();}


FittedParameter::FittedParameter(std::string name, double val, Limit error) : name(name), val(val), error(error) {}
FittedParameter::FittedParameter(std::string name, double val, double error) : name(name), val(val), error({-error, +error}) {}
FittedParameter::FittedParameter(const Parameter& param, double val, Limit error) : name(param.name), val(val), error(error) {}
FittedParameter::FittedParameter(const Parameter& param, double val, double error) : name(param.name), val(val), error(-error, +error) {}