// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::rigidbody::parameter::decay {
    /**
     * @brief The strategy for decaying the amplitude of a parameter.
     */
    class DecayStrategy {
        public:
            DecayStrategy() = default;
            DecayStrategy(unsigned int iterations) : iterations(iterations) {}
            virtual ~DecayStrategy() = default;

            /**
             * @brief Get the factor by which to decay the parameter and increment the internal counter. 
             */
            virtual double next() = 0;

            void set_iterations(unsigned int iterations) {
                this->iterations = iterations;
                set_characteristic_time(iterations);
            }

            int get_iterations() const {return iterations;}

            /**
             * @brief Restart the decay from its initial value.
             */
            void reset() {draws = 0;}

        protected:
            /**
             * @brief Consume a single draw and return its index.
             */
            unsigned int next_draw();

            unsigned int iterations = 0;
            unsigned int draws = 0;

        private:
            bool overdrawn_warning_issued = false;

            /**
             * @brief Set the characteristic time scale for this decay strategy.
             */
            virtual void set_characteristic_time(unsigned int iterations) = 0;
    };
}