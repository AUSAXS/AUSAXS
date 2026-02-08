// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <settings/RigidBodySettings.h>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<rigidbody::transform::TransformStrategy> create_transform_strategy(observer_ptr<Rigidbody> body);
    std::unique_ptr<rigidbody::transform::TransformStrategy> create_transform_strategy(observer_ptr<Rigidbody> body, settings::rigidbody::TransformationStrategyChoice choice);
}