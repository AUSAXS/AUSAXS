// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>

namespace ausaxs::rigidbody::parameter::decay {
    /**
     * @brief A DecayStrategy that does not decay the amplitudes at all. 
     */
    class NoDecay : public DecayStrategy {
        public:
            NoDecay() = default;
            ~NoDecay() = default;

            double next() override {return 1;}

        private:
            void set_characteristic_time(unsigned int) override {}
    };
}