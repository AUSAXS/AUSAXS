// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/decay/DecayStrategy.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> create_decay_strategy(int iterations);
    std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> create_decay_strategy(int iterations, settings::rigidbody::DecayStrategyChoice choice);
}