/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/SimpleParameterGeneration.h>
#include <rigidbody/parameters/TranslationsOnly.h>
#include <rigidbody/parameters/RotationsOnly.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace rigidbody::parameter;

std::unique_ptr<ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(unsigned int iterations, double translate_amp, double rotate_amp) {
    return create_parameter_strategy(rigidbody::factory::create_decay_strategy(iterations), translate_amp, rotate_amp, settings::rigidbody::parameter_generation_strategy);
}

std::unique_ptr<ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(unsigned int iterations, double translate_amp, double rotate_amp, settings::rigidbody::ParameterGenerationStrategyChoice choice) {
    return create_parameter_strategy(rigidbody::factory::create_decay_strategy(iterations), translate_amp, rotate_amp, choice);
}

std::unique_ptr<ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> decay_strategy, double translate_amp, double rotate_amp, settings::rigidbody::ParameterGenerationStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::ParameterGenerationStrategyChoice::Simple:
            return std::make_unique<SimpleParameterGeneration>(std::move(decay_strategy), translate_amp, rotate_amp);
        case settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly:
            return std::make_unique<RotationsOnly>(std::move(decay_strategy), translate_amp, rotate_amp);
        case settings::rigidbody::ParameterGenerationStrategyChoice::TranslationsOnly:
            return std::make_unique<TranslationsOnly>(std::move(decay_strategy), translate_amp, rotate_amp);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_parameter_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}