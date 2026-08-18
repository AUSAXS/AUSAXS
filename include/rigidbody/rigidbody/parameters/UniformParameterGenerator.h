// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief Draws every enabled parameter component from a uniform distribution spanning its amplitude.
     */
    class UniformParameterGenerator : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~UniformParameterGenerator() override = default;

            BodyTransformParametersRelative next(int ibody) override;
    };
}
