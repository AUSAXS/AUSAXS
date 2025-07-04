// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/parameters/decay/ExponentialDecay.h>
#include <rigidbody/parameters/decay/LinearDecay.h>
#include <rigidbody/parameters/decay/NoDecay.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> rigidbody::factory::create_decay_strategy(unsigned int iterations) {
    return create_decay_strategy(iterations, settings::rigidbody::decay_strategy);
}

std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> rigidbody::factory::create_decay_strategy(unsigned int iterations, settings::rigidbody::DecayStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::DecayStrategyChoice::Linear:
            return std::make_unique<rigidbody::parameter::decay::LinearDecay>(iterations);
        case settings::rigidbody::DecayStrategyChoice::Exponential:
            return std::make_unique<rigidbody::parameter::decay::ExponentialDecay>(iterations);
        case settings::rigidbody::DecayStrategyChoice::None:
            return std::make_unique<rigidbody::parameter::decay::NoDecay>();
        default:
            throw except::unknown_argument("rigidbody::parameters::factory::create_decay_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}