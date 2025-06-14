// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::selection;

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(observer_ptr<const Rigidbody> body) {
    return create_selection_strategy(body, settings::rigidbody::body_select_strategy);
}

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(observer_ptr<const Rigidbody> body, settings::rigidbody::BodySelectStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::BodySelectStrategyChoice::RandomBodySelect:
            return std::make_unique<RandomBodySelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect:
            return std::make_unique<RandomConstraintSelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::SequentialBodySelect:
            return std::make_unique<SequentialBodySelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::SequentialConstraintSelect:
            return std::make_unique<SequentialConstraintSelect>(body);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_selection_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}