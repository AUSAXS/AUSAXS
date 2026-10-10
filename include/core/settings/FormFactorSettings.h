// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>

namespace ausaxs::settings {
    /// @brief Settings controlling which form factor types get a slot of their own.
    struct EXPORT form_factor {
        // The maximum number of active form factor slots, including EXCLUDED_VOLUME, WATER and OTHER, which are always present.
        // The partial histograms scale quadratically with this number, so it is the dominant memory cost of the partial managers.
        // Excess types are folded onto OTHER in order of increasing abundance.
        static int max_types;

        // The minimum fraction of the atoms a form factor type must make up to get a slot of its own. Rarer types are folded onto OTHER.
        static double min_fraction;
    };
}
