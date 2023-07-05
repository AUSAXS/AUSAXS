#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>

#include <memory>

namespace settings::rigidbody {enum class ParameterGenerationStrategyChoice;}
namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp);

        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp, const settings::rigidbody::ParameterGenerationStrategyChoice& choice);
    }
}