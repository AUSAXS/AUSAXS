#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>

namespace ausaxs::rigidbody::parameter::decay {
    class ExponentialDecay : public DecayStrategy {
        public:
            ExponentialDecay(unsigned int max_iterations);
            ~ExponentialDecay();

            double next() override;

            /**
             * @brief Set the number of iterations required to reach an amplitude of 1/e.
             */
            void set_characteristic_time(unsigned int iterations) override;

        private:
            double decay_rate;
    };
}