#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/SimpleParameterGeneration.h>
#include <rigidbody/parameters/RotationsOnly.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

std::unique_ptr<rigidbody::ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(unsigned int iterations, double translate_amp, double rotate_amp) {
    return create_parameter_strategy(iterations, translate_amp, rotate_amp, settings::rigidbody::parameter_generation_strategy);
}

std::unique_ptr<rigidbody::ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(unsigned int iterations, double translate_amp, double rotate_amp, const settings::rigidbody::ParameterGenerationStrategyChoice& choice) {
    switch (choice) {
        case settings::rigidbody::ParameterGenerationStrategyChoice::Simple:
            return std::make_unique<rigidbody::SimpleParameterGeneration>(iterations, translate_amp, rotate_amp);
        case settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly:
            return std::make_unique<rigidbody::RotationsOnly>(iterations, translate_amp, rotate_amp);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_parameter_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}