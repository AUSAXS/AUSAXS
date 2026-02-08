// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/RigidBodySettings.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> create_decay_strategy(unsigned int iterations);
    std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> create_decay_strategy(unsigned int iterations, settings::rigidbody::DecayStrategyChoice choice);
}