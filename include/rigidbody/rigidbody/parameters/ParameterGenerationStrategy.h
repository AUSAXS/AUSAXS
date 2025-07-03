// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/MathFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/parameters/Parameter.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>
#include <utility/observer_ptr.h>

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
            ParameterGenerationStrategy(observer_ptr<const RigidBody> molecule, unsigned int iterations, double length_start, double rad_start);

            /**
             * @brief Construct a new parameter generation strategy.
             * 
             * @param decay_strategy The decay strategy to use.
             * @param length_start The start length of the generated translation vectors. 
             * @param rad_start The start angle in radians of the generated rotations. 
             */
            ParameterGenerationStrategy(observer_ptr<const RigidBody> molecule, std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy, double length_start, double rad_start);

            virtual ~ParameterGenerationStrategy();

            /**
             * @brief Get the next parameter in the sequence for the given body.
             */
            virtual Parameter next(int ibody) = 0;

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
            observer_ptr<const RigidBody> molecule;
            std::mt19937 generator;
            std::uniform_real_distribution<double> translation_dist; // Random number distribution for translations. 
            std::uniform_real_distribution<double> rotation_dist;    // Random number distribution for rotations. 
            std::uniform_real_distribution<double> symmetry_dist;    // Random number distribution for symmetry transforms.
            std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy;
    };
}