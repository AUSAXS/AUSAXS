// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/Exceptions.h>
#include <utility/Console.h>

using namespace ausaxs::except;

base::base(const char* msg) : msg(msg) {
    console::print_critical(msg);
}

base::base(const std::string msg) : msg(msg) {
    console::print_critical(msg);
}

const char* base::what() const noexcept {return msg.data();}