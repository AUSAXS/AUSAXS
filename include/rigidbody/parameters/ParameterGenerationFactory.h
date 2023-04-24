#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp, settings::rigidbody::ParameterGenerationStrategyChoice choice = settings::rigidbody::parameter_generation_strategy);
    }
}