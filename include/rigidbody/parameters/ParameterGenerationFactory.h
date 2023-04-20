#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <utility/SmartOption.h>

namespace settings {
    namespace rigidbody {
        enum class ParameterGenerationStrategyChoice {
            Simple,         // Generate translation and rotation parameters. Their amplitudes decays linearly with the iteration number.
            RotationsOnly   // Only generate rotation parameters. The amplitudes decays linearly with the iteration number.
        };

        extern settings::detail::SmartOption<ParameterGenerationStrategyChoice> parameter_generation_strategy;
    }
}

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotation_amp, settings::rigidbody::ParameterGenerationStrategyChoice choice = settings::rigidbody::parameter_generation_strategy);
    }
}