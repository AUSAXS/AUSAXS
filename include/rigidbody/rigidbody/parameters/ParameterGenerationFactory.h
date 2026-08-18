// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/RigidbodyFwd.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    /**
     * @brief Zero the amplitudes of every component the given choice excludes, leaving the rest untouched.
     */
    parameter::ParameterAmplitudes restrict_to(parameter::ParameterAmplitudes amplitudes, settings::rigidbody::ParameterGenerationStrategyChoice choice);

    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(
        observer_ptr<const Rigidbody> molecule, unsigned int iterations, const parameter::ParameterAmplitudes& amplitudes
    );

    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(
        observer_ptr<const Rigidbody> molecule, std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> decay_strategy,
        const parameter::ParameterAmplitudes& amplitudes
    );

    /**
     * @brief Create a parameter strategy with the components excluded by the given choice disabled.
     *
     * Used for the default, unscripted strategy; a script instead resolves its "mode" argument into the amplitudes directly.
     */
    std::unique_ptr<parameter::ParameterGenerationStrategy> create_parameter_strategy(
        observer_ptr<const Rigidbody> molecule, unsigned int iterations, const parameter::ParameterAmplitudes& amplitudes,
        settings::rigidbody::ParameterGenerationStrategyChoice choice
    );
}
