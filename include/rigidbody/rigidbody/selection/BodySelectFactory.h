// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace ausaxs::rigidbody {
    namespace factory {
        /**
         * @brief Prepare a body selection strategy.
         */
        std::unique_ptr<selection::BodySelectStrategy> create_selection_strategy(observer_ptr<const RigidBody> molecule);

        /**
         * @brief Prepare a body selection strategy.
         */
        std::unique_ptr<selection::BodySelectStrategy> create_selection_strategy(observer_ptr<const RigidBody> molecule, settings::rigidbody::BodySelectStrategyChoice choice);
    }
}