// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

ParameterAmplitudes rigidbody::factory::restrict_to(ParameterAmplitudes amplitudes, settings::rigidbody::ParameterGenerationStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::ParameterGenerationStrategyChoice::Simple:
            return amplitudes;
        case settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly:
            return {.rotation = amplitudes.rotation};
        case settings::rigidbody::ParameterGenerationStrategyChoice::TranslationsOnly:
            return {.translation = amplitudes.translation};
        case settings::rigidbody::ParameterGenerationStrategyChoice::SymmetryOnly:
            return {.symmetry_translation = amplitudes.symmetry_translation, .symmetry_rotation = amplitudes.symmetry_rotation};
        default: 
            throw except::unknown_argument("rigidbody::factory::restrict_to: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}

std::unique_ptr<ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(
    observer_ptr<const Rigidbody> molecule, unsigned int iterations, const ParameterAmplitudes& amplitudes
) {
    return create_parameter_strategy(molecule, rigidbody::factory::create_decay_strategy(iterations), amplitudes);
}

std::unique_ptr<ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(
    observer_ptr<const Rigidbody> molecule, std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> decay_strategy, const ParameterAmplitudes& amplitudes
) {
    return std::make_unique<UniformParameterGenerator>(molecule, std::move(decay_strategy), amplitudes);
}

std::unique_ptr<ParameterGenerationStrategy> rigidbody::factory::create_parameter_strategy(
    observer_ptr<const Rigidbody> molecule, unsigned int iterations, settings::rigidbody::ParameterGenerationStrategyChoice choice
) {
    return create_parameter_strategy(molecule, rigidbody::factory::create_decay_strategy(iterations), restrict_to(default_amplitudes(molecule), choice));
}
