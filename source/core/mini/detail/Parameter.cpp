// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/detail/Parameter.h>
#include <mini/detail/FittedParameter.h>

using namespace ausaxs::mini;

Parameter::Parameter() = default;

Parameter::Parameter(const std::string& name, const Limit& bounds) noexcept: name(name), bounds(bounds) {}

Parameter::Parameter(const std::string& name, double guess) noexcept: Parameter(name, guess, {0, 0}) {}

Parameter::Parameter(const std::string& name, double guess, const Limit& bounds) noexcept: name(name), guess(guess), bounds(bounds) {}

Parameter::Parameter(const mini::FittedParameter& p) noexcept {
    *this = p;
}

bool Parameter::has_bounds() const noexcept {
    return bounds.has_value();
}

bool Parameter::has_guess() const noexcept {
    return guess.has_value();
}

bool Parameter::has_name() const noexcept {
    return !name.empty();
}

bool Parameter::empty() const noexcept {
    return !(has_name() && (has_bounds() || has_guess()));
}

std::string Parameter::to_string() const {
    std::string s = name;
    if (has_guess()) {s += " guess " + std::to_string(guess.value());}
    if (has_bounds()) {s += " bounds [" + std::to_string(bounds.value().min) + std::to_string(bounds.value().max) + "]";}
    return s;
}

Parameter& Parameter::operator=(const mini::FittedParameter& other) noexcept {
    name = other.name;
    guess = other.value;
    bounds = other.error + other.value;
    return *this;
}