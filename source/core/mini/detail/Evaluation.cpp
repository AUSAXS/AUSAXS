// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/detail/Evaluation.h>

#include <string>

using namespace ausaxs;

std::string mini::Evaluation::to_string() const {
    std::string s;
    for (auto val : vals) {
        s += std::to_string(val) + " ";
    }
    s += std::to_string(fval);
    return s;
}