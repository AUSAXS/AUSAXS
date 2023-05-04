#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>

namespace rigidbody::parameters::decay {
    class LinearDecay : public DecayStrategy {
        public:
            LinearDecay(unsigned int max_iterations);
            LinearDecay(double decay_rate);
            ~LinearDecay() = default;

            double get_factor() override;

            /**
             * @brief Set the number of iterations required to reach an amplitude of 0.5.
             */
            void set_characteristic_time(unsigned int iterations);
        
        private:
            double decay_rate = 0.01;
    };
}