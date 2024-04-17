/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace rigidbody::selection;

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(const rigidbody::RigidBody* body) {
    return create_selection_strategy(body, settings::rigidbody::body_select_strategy);
}

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(const rigidbody::RigidBody* body, settings::rigidbody::BodySelectStrategyChoice choice) {
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