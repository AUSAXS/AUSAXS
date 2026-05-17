// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::settings {
    namespace form_factor {
        /**
        * @brief The compile-time upper bound on the number of active form factor types (including EXV).
        *        Increasing this value allows more atom-type classes to be supported at the cost of slightly larger histograms.
        */
        constexpr int max_ff_types = 15;
    };
}