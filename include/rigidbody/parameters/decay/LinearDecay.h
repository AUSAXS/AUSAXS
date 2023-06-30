#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>

namespace rigidbody::parameters::decay {
    /**
     * @brief A DecayStrategy that decays the factor linearly.
     */
    class LinearDecay : public DecayStrategy {
        public:
            /**
             * @brief Create a new LinearDecay object.
             *        The factor will be zero after the given number of maximum iterations.
             */
            LinearDecay(unsigned int max_iterations);
            ~LinearDecay();

            double get_factor() override;
        
        private:
            /**
             * @brief Set the number of iterations required to reach an amplitude of 0.5.
             */
            void set_characteristic_time(unsigned int iterations);

            double decay_rate = 0;
    };
}