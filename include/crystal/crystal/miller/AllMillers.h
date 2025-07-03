// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>

namespace ausaxs::crystal {
    /**
     * @brief Generates all miller indices within the range specified by
     *          settings::crystal::h
     *          settings::crystal::k
     *          settings::crystal::l
     * 
     *        The maximum allowed length of the indices is specified by
     *          settings::crystal::max_q
     *        which guarantees that the indices spans a spherical volume.
     */
    class AllMillers : public MillerGenerationStrategy {
        public: 
            AllMillers(unsigned int h, unsigned int k, unsigned int l);
            ~AllMillers() override = default;

            std::vector<Miller> generate() const override;
        private: 
            int h, k, l;
    };
}