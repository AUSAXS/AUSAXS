// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>

namespace ausaxs::crystal {

    /**
     * @brief Generates a subset of the indices within the range specified by
     *           settings::crystal::h
     *           settings::crystal::k
     *           settings::crystal::l
     *        All independent basis vectors out to the length 
     *           settings::crystal::reduced::basis_q
     *        are generated. These basis vectors are then multiplied by a constant
     *        to generate the full set of indices within the range specified by
     *           settings::crystal::max_q
     *        This is significantly faster than using the AllMillers strategy, but also loses some information.
     */
    class ReducedMillers : public MillerGenerationStrategy {
        public:
            ReducedMillers();
            ReducedMillers(unsigned int h, unsigned int k, unsigned int l);
            ~ReducedMillers() override = default;

            std::vector<Miller> generate() const override;

        protected:
            std::vector<Miller> generate_independent_bases(double limit = -1) const;

        private:
            int h, k, l;
    };
}