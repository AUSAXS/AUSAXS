/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <mini/detail/Parameter.h>
#include <mini/detail/FittedParameter.h>

using namespace mini;

Parameter::Parameter() noexcept = default;

Parameter::Parameter(const std::string& name, const Limit& bounds) noexcept: name(name), bounds(bounds) {}

Parameter::Parameter(const std::string& name, double guess) noexcept: Parameter(name, guess, {0, 0}) {}

Parameter::Parameter(const std::string& name, double guess, const Limit& bounds) noexcept: name(name), guess(guess), bounds(bounds) {}

Parameter::Parameter(const mini::FittedParameter& p) noexcept {
    *this = p;
}

Parameter::~Parameter() = default;

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

std::string Parameter::to_string() const noexcept {
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