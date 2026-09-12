// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <math/Exceptions.h>

#include <utility/Console.h>

using namespace ausaxs::math::except;

base::base(const char* msg) : msg(msg) {
    console::print_critical(msg);
}

base::base(std::string msg) : msg(std::move(msg)) {
    console::print_critical(this->msg);
}

const char* base::what() const noexcept {return msg.data();}