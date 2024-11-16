#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace ausaxs::rigidbody::parameter {
    /**
     * @brief Thread-safe parameter generation strategy. The current step size scales linearly with the iteration number. 
     */
    class SimpleParameterGeneration : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~SimpleParameterGeneration() override;

            Parameter next() override;
    };
}