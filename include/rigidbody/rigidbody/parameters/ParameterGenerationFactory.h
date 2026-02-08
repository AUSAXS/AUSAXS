// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/RigidbodyFwd.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(
        observer_ptr<const Rigidbody> molecule, unsigned int iterations, double translate_amp, double rotation_amp
    );

    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(
        observer_ptr<const Rigidbody> molecule, unsigned int iterations, double translate_amp, double rotation_amp, 
        settings::rigidbody::ParameterGenerationStrategyChoice choice
    );

    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(
        observer_ptr<const Rigidbody> molecule, std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> decay_strategy, 
        double translate_amp, double rotate_amp, settings::rigidbody::ParameterGenerationStrategyChoice choice
    );
}