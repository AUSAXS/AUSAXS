#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

namespace rigidbody::parameter {
    class TranslationsOnly : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~TranslationsOnly() override;

            Parameter next() override;
    };
}