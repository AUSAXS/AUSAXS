// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>

namespace ausaxs::rigidbody::parameter::decay {
    class ExponentialDecay : public DecayStrategy {
        public:
            ExponentialDecay(int max_iterations);
            ~ExponentialDecay() override;

            double next() override;

            /**
             * @brief Set the number of iterations required to reach an amplitude of 1/e.
             */
            void set_characteristic_time(int iterations) override;

        private:
            double decay_rate;
    };
}