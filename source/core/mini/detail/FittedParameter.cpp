// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/detail/FittedParameter.h>
#include <mini/detail/Parameter.h>

using namespace ausaxs::mini;

FittedParameter::FittedParameter(const std::string& name, double val, const Limit& error) noexcept : name(name), value(val), error(error) {}

FittedParameter::FittedParameter(const std::string& name, double val, double error) noexcept: name(name), value(val), error({-error, +error}) {}

FittedParameter::FittedParameter(const Parameter& param, double val, const Limit& error) noexcept: name(param.name), value(val), error(error) {}

FittedParameter::FittedParameter(const Parameter& param, double val, double error) noexcept: name(param.name), value(val), error(-error, +error) {}

std::string FittedParameter::to_string() const noexcept {
    return name + " " + std::to_string(value) + " " + error.to_string();
}

double FittedParameter::mean_error() const noexcept {
    return (value - error.min)*0.5 + (error.max - value)*0.5;
}

void FittedParameter::set_error(double error) noexcept {
    this->error = Limit(-error, error);
}

void FittedParameter::set_error(double min, double max) noexcept {
    this->error = Limit(min, max);
}