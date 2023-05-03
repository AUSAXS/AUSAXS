#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomSelect.h>
#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <utility/Exceptions.h>

std::unique_ptr<rigidbody::selection::BodySelectStrategy> rigidbody::factory::create_selection_strategy(const rigidbody::RigidBody* body, settings::rigidbody::BodySelectStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::BodySelectStrategyChoice::RandomSelect:
            return std::make_unique<rigidbody::selection::RandomSelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect:
            return std::make_unique<rigidbody::selection::RandomConstraintSelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::SequentialSelect:
            return std::make_unique<rigidbody::selection::SequentialSelect>(body);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_selection_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}