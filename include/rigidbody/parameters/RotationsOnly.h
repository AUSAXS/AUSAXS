#pragma once

#include <utility/Concepts.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace rigidbody::parameter {
    /**
     * @brief Thread-safe parameter generation strategy. The current step size scales linearly with the iteration number. 
     */
    class RotationsOnly : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~RotationsOnly() override;

            Parameter next() override;
    };
}