// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <mini/detail/Landscape.h>

namespace ausaxs::em {
    namespace detail {
        /**
         * @brief Specialized landscape class for the ImageStack class. 
         * 
         * Since the fit is performed as two 1D fits, we only get a single "strip" of the landscape for each iteration. 
         * This class extends the mini::Landscape class to store both the strip and the cutoff value corresponding to that iteration.
         */
        struct ExtendedLandscape {
            /**
             * @brief Constructor.
             */
            ExtendedLandscape(double cutoff, double mass, double volume, mini::Landscape&& l) : cutoff(cutoff), mass(mass), volume(volume), strip(std::move(l)) {}

            double cutoff, mass, volume;
            mini::Landscape strip;
        };
    }
}