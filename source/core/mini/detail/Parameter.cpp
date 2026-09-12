// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/detail/Parameter.h>

#include <mini/detail/FittedParameter.h>

using namespace ausaxs::mini;

Parameter::Parameter() = default;

Parameter::Parameter(std::string  name, const Limit& bounds) noexcept: name(std::move(name)), bounds(bounds) {}

Parameter::Parameter(std::string name, double guess) noexcept: Parameter(std::move(name), guess, {0, 0}) {}

Parameter::Parameter(std::string  name, double guess, const Limit& bounds) noexcept: name(std::move(name)), guess(guess), bounds(bounds) {}

Parameter::Parameter(const mini::FittedParameter& p) noexcept {
    *this = p;
}

bool Parameter::empty() const noexcept {
    return !(has_name() && (has_bounds() || has_guess()));
}

std::string Parameter::to_string() const {
    std::string s = name;
    if (guess.has_value()) {s += " guess " + std::to_string(*guess);}
    if (bounds.has_value()) {s += " bounds [" + std::to_string(bounds->min) + std::to_string(bounds->max) + "]";}
    return s;
}

Parameter& Parameter::operator=(const mini::FittedParameter& other) noexcept {
    name = other.name;
    guess = other.value;
    bounds = other.error + other.value;
    return *this;
}