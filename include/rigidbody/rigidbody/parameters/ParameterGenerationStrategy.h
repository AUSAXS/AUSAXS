#pragma once

#include <math/MathFwd.h>
#include <utility/Concepts.h>
#include <rigidbody/parameters/Parameter.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>

#include <random>

namespace ausaxs::rigidbody::parameter {    
    /**
     * @brief Thread-safe superclass for parameter generation strategies.
     */
    class ParameterGenerationStrategy {
        public:
            /**
             * @brief Construct a new parameter generation strategy.
             * 
             * @param iterations The expected number of iterations. This is used to determine the linear decay rate.
             * @param length_start The start length of the generated translation vectors. 
             * @param rad_start The start angle in radians of the generated rotations. 
             */
            ParameterGenerationStrategy(unsigned int iterations, double length_start, double rad_start);

            /**
             * @brief Construct a new parameter generation strategy.
             * 
             * @param decay_strategy The decay strategy to use.
             * @param length_start The start length of the generated translation vectors. 
             * @param rad_start The start angle in radians of the generated rotations. 
             */
            ParameterGenerationStrategy(std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy, double length_start, double rad_start);

            /**
             * @brief Destructor.
             */
            virtual ~ParameterGenerationStrategy();

            /**
             * @brief Get the next parameter in the sequence. 
             */
            virtual Parameter next() = 0;

            /**
             * @brief Set the decay strategy.
             */
            void set_decay_strategy(std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy);

            /**
             * @brief Set the maximum translation distance to the given value.
             */
            void set_max_translation_distance(double distance);

            /**
             * @brief Set the maximum rotation angle to the given value.
             */
            void set_max_rotation_angle(double radians);

        protected:
            std::mt19937 generator;                                             // The random number generator. 
            std::uniform_real_distribution<double> translation_dist;            // The random number distribution for translations. 
            std::uniform_real_distribution<double> rotation_dist;               // The random number distribution for rotations. 
            std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy;    // The decay strategy.
    };
}