#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp);

        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp, const settings::rigidbody::ParameterGenerationStrategyChoice& choice);
    }
}