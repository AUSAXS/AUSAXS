// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <array>
#include <numbers>

namespace ausaxs::constants::form_factor {
    namespace {
        constexpr double s_to_q_factor = 1./(4*4*std::numbers::pi*std::numbers::pi); // q = 4πs --> s = q/(4pi)

        /**
         * @brief Convert a scattering vector from s to q.
         *        This is purely for convenience, such that the tabulated values are easier to read.
         */
        constexpr std::array<double, 5> s_to_q(std::array<double, 5> a) {
            for (unsigned int i = 0; i < 5; ++i) {
                a[i] *= s_to_q_factor;
            }
            return a;
        }
    }
}