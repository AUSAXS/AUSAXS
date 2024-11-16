#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    /**
     * @brief Prepare a decay class.
     */
    std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> create_decay_strategy(unsigned int iterations);

    /**
     * @brief Prepare a decay class.
     */
    std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> create_decay_strategy(unsigned int iterations, settings::rigidbody::DecayStrategyChoice choice);
}