// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <shell/Option.h>

using namespace ausaxs::shell;

Option::Option() noexcept = default;
Option::Option(const std::string& name, const std::string& value) : name(name), value(value) {}
Option::~Option() = default;

std::string Option::get() const {
    if (value.empty()) {
        return name;
    }
    return name + " " + value;
}
Argument::Argument(const std::string& name, const std::string& value) : Option(name, value) {}
Argument::Argument(const std::string& name, double value) : Option(name, std::to_string(value)) {}
Argument::Argument(const std::string& name, int value) : Option(name, std::to_string(value)) {}
Argument::Argument(const std::string& name, unsigned int value) : Option(name, std::to_string(value)) {}

Flag::Flag(const std::string& name) : Option(name, "") {}