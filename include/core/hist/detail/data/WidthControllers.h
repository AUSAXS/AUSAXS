// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <settings/Flags.h>

namespace ausaxs::hist::detail {
    struct ConstantWidth {
        static consteval float get() {return 1./ausaxs::constants::axes::d_axis.width();}
    };

    struct VariableWidth {
        static float get() {return settings::flags::inv_bin_width;}
    };

    template<bool variable_bin_width>
    struct WidthController {
        constexpr static float get_inv_width() {
            if constexpr (variable_bin_width) {
                return VariableWidth::get();
            } else {
                return ConstantWidth::get();
            }
        }
    };
}