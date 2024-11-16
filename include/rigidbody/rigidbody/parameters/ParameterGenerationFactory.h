#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    /**
     * @brief Prepare a constraint generator. 
     */
    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp);

    /**
     * @brief Prepare a constraint generator. 
     */
    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp, settings::rigidbody::ParameterGenerationStrategyChoice choice);

    /**
     * @brief Prepare a constraint generator. 
     */
    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> decay_strategy, double translate_amp, double rotate_amp, settings::rigidbody::ParameterGenerationStrategyChoice choice);
}