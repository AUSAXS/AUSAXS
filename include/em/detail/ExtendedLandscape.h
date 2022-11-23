#pragma once

#include <mini/utility/Landscape.h>

namespace em {
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
            ExtendedLandscape(double cutoff, mini::Landscape&& l) : cutoff(cutoff), strip(std::move(l)) {}

            /**
             * @brief Constructor.
             * 
             * Takes ownership of the given landscape.
             */
            ExtendedLandscape(double cutoff, mini::Landscape& l) : cutoff(cutoff), strip(std::move(l)) {}

            double cutoff;
            mini::Landscape strip;
        };
    }
}