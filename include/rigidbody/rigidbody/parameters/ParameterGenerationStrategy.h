// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/MathFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/parameters/ParameterAmplitudes.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>
#include <utility/observer_ptr.h>

#include <random>

namespace ausaxs::rigidbody::parameter {    
    /**
     * @brief Superclass for parameter generation strategies.
     * 
     * Generates relative transformations to be applied to bodies or rigid groups.
     */
    class ParameterGenerationStrategy {
        public:
            /**
             * @brief Construct a new parameter generation strategy.
             * 
             * @param iterations The expected number of iterations. This is used to determine the linear decay rate.
             * @param amplitudes The maximum amplitude of each parameter component. Components with a zero amplitude are not generated.
             */
            ParameterGenerationStrategy(observer_ptr<const Rigidbody> molecule, unsigned int iterations, const ParameterAmplitudes& amplitudes);

            /**
             * @brief Construct a new parameter generation strategy.
             * 
             * @param decay_strategy The decay strategy to use.
             * @param amplitudes The maximum amplitude of each parameter component. Components with a zero amplitude are not generated.
             */
            ParameterGenerationStrategy(observer_ptr<const Rigidbody> molecule, std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy, const ParameterAmplitudes& amplitudes);

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
             * @brief Set the maximum translation distance to the given value. Zero disables body translations.
             */
            void set_max_translation_distance(double distance);

            /**
             * @brief Set the maximum rotation angle to the given value. Zero disables body rotations.
             */
            void set_max_rotation_angle(double radians);

            /**
             * @brief Set the maximum symmetry offset translation distance to the given value. Zero disables symmetry translations.
             */
            void set_max_symmetry_translation_distance(double distance);

            /**
             * @brief Set the maximum symmetry frame reorientation to the given value. Zero disables symmetry rotations.
             */
            void set_max_symmetry_rotation_angle(double radians);

            const ParameterAmplitudes& get_amplitudes() const;

        protected:
            // Draw a uniform deviate in [-amplitude, amplitude].
            double draw(double amplitude);

            // Draw a unit vector pointing in a uniformly distributed direction.
            Vector3<double> draw_direction();

            observer_ptr<const Rigidbody> rigidbody;
            ParameterAmplitudes amplitudes;
            std::unique_ptr<parameter::decay::DecayStrategy> decay_strategy;

        private:
            // All components are drawn from this single distribution and scaled by their amplitude.
            std::uniform_real_distribution<double> unit_dist{-1, 1};
    };
}
