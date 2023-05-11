#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>

#include <memory>

namespace settings::rigidbody {enum class DecayStrategyChoice;}
namespace rigidbody::parameters::factory {
    /**
     * @brief Prepare a decay class.
     */
    std::unique_ptr<rigidbody::parameters::decay::DecayStrategy> create_decay_strategy(unsigned int iterations);

    /**
     * @brief Prepare a decay class.
     */
    std::unique_ptr<rigidbody::parameters::decay::DecayStrategy> create_decay_strategy(unsigned int iterations, const settings::rigidbody::DecayStrategyChoice& choice);
}