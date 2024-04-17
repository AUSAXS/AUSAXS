#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace rigidbody::parameter {
    class RotationsOnly : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~RotationsOnly() override;

            Parameter next() override;
    };
}