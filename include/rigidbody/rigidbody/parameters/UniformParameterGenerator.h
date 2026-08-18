// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief Draws every enabled parameter component from a uniform distribution spanning its amplitude.
     *
     * Which components are generated follows from the amplitudes alone: a component with a zero amplitude is
     * left unset rather than being generated as a zero delta.
     */
    class UniformParameterGenerator : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~UniformParameterGenerator() override = default;

            BodyTransformParametersRelative next(int ibody) override;
    };
}
