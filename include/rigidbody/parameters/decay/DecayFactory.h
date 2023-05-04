#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace rigidbody::parameters::factory {
    /**
     * @brief Prepare a decay class.
     */
    std::unique_ptr<rigidbody::parameters::decay::DecayStrategy> create_decay_strategy(settings::rigidbody::DecayStrategyChoice choice = settings::rigidbody::decay_strategy);
}