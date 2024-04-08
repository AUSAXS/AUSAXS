/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomSelect.h>
#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace rigidbody::selection;

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(const rigidbody::RigidBody* body) {
    return create_selection_strategy(body, settings::rigidbody::body_select_strategy);
}

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(const rigidbody::RigidBody* body, settings::rigidbody::BodySelectStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::BodySelectStrategyChoice::RandomSelect:
            return std::make_unique<RandomSelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect:
            return std::make_unique<RandomConstraintSelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::SequentialSelect:
            return std::make_unique<SequentialSelect>(body);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_selection_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}