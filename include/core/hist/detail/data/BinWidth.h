// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/InternalState.h>

namespace ausaxs::hist::detail {
    /**
     * @brief The inverse width of the distance bins, at the precision the distance kernels bin with.
     *
     * Read it once, outside any loop that writes a histogram, and pass the value on: the setting is a global double, so the
     * compiler must assume that every store into a double histogram may overwrite it, and would read it again after each one.
     */
    inline float inv_bin_width() {return static_cast<float>(settings::internal_state::inv_bin_width);}
}
