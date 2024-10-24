/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <mini/detail/Result.h>
#include <utility/Exceptions.h>
#include <mini/detail/FittedParameter.h>

#include <algorithm>
#include <utility>

using namespace ausaxs::mini;

Result::Result(const FittedParameter& param, double fval, unsigned int fevals) noexcept : parameters({param}), fval(fval), fevals(fevals) {}

Result::Result(const std::vector<FittedParameter>& params, double fval, unsigned int fevals) noexcept : parameters(params), fval(fval), fevals(fevals) {}

void Result::add_parameter(const FittedParameter& param) noexcept {
    parameters.push_back(param);
}

unsigned int Result::size() const noexcept {
    return parameters.size();
}

unsigned int Result::dim() const noexcept {
    return size();
}

const FittedParameter& Result::get_parameter(const std::string& name) const {
    auto pos = std::find_if(parameters.begin(), parameters.end(), [&name] (const FittedParameter& param) {return param.name == name;});
    if (pos == parameters.end()) {throw except::unknown_argument("Result::get_parameter: No parameter named \"" + name + "\" was found.");}
    return *pos;
}

FittedParameter& Result::get_parameter(const std::string& name) {
    return const_cast<FittedParameter&>(std::as_const(*this).get_parameter(name));
}

const FittedParameter& Result::get_parameter(unsigned int index) const {
    if (size() < index) {throw except::out_of_bounds("Result::get_parameter: Index \"" + std::to_string(index) + "\" is out of bounds (" + std::to_string(size()) + ").");}
    return parameters[index];
}

FittedParameter& Result::get_parameter(unsigned int index) {
    return const_cast<FittedParameter&>(std::as_const(*this).get_parameter(index));
}

const std::vector<FittedParameter>& Result::get_parameters() const {
    return parameters;
}

std::vector<FittedParameter> Result::get_parameters() {
    return parameters;
}

std::vector<double> Result::get_parameter_values() const {
    std::vector<double> values;
    for (const auto& p : parameters) {
        values.push_back(p.value);
    }
    return values;
}

const FittedParameter& Result::operator[](unsigned int index) const {
    return get_parameter(index);
}

FittedParameter& Result::operator[](unsigned int index) {
    return get_parameter(index);
}