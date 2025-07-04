// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <settings/RigidBodySettings.h>

namespace ausaxs::rigidbody::factory {
    /**
     * @brief Prepare a transformation strategy.
     */
    std::unique_ptr<rigidbody::transform::TransformStrategy> create_transform_strategy(RigidBody* body);

    /**
     * @brief Prepare a transformation strategy.
     */
    std::unique_ptr<rigidbody::transform::TransformStrategy> create_transform_strategy(RigidBody* body, const settings::rigidbody::TransformationStrategyChoice& choice);
}