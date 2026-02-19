// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/MathFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>
#include <utility/observer_ptr.h>

#include <random>

namespace ausaxs::rigidbody::parameter {    
    /**
     * @brief Thread-safe superclass for parameter generation strategies.
     * 
     * Generates relative transformations to be applied to bodies or rigid groups.
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
            ParameterGenerationStrategy(observer_ptr<const Rigidbody> molecule, unsigned int iterations, double length_start, double rad_start);

            /**
             * @brief Construct a new parameter generation strategy.
             * 
             * @param decay_strategy The decay strategy to use.
             * @param length_start The start length of the generated translation vectors. 
             * @param rad_start The start angle in radians of the generated rotations. 
             */
            ParameterGenerationStrategy(observer_ptr<const Rigidbody> molecule, std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy, double length_start, double rad_start);

            virtual ~ParameterGenerationStrategy();

            /**
             * @brief Get the next relative transformation for the given body.
             */
            virtual BodyTransformParametersRelative next(int ibody) = 0;

            /**
             * @brief Set the decay strategy.
             */
            void set_decay_strategy(std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy);
            observer_ptr<parameter::decay::DecayStrategy> get_decay_strategy() const;

            /**
             * @brief Set the maximum translation distance to the given value.
             */
            void set_max_translation_distance(double distance);

            /**
             * @brief Set the maximum rotation angle to the given value.
             */
            void set_max_rotation_angle(double radians);

        protected:
            observer_ptr<const Rigidbody> rigidbody;
            std::uniform_real_distribution<double> translation_dist;
            std::uniform_real_distribution<double> rotation_dist;
            std::uniform_real_distribution<double> translation_symmetry_dist;
            std::uniform_real_distribution<double> rotation_symmetry_dist;
            std::uniform_real_distribution<double> angle_symmetry_dist;
            std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy;
    };
}