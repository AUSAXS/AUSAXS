// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <rigidbody/selection/ParameterMaskStrategy.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::selection;

namespace {
    std::unique_ptr<ParameterMaskStrategy> create_mask_strategy(settings::rigidbody::ParameterMaskStrategyChoice choice) {
        switch (choice) {
            case settings::rigidbody::ParameterMaskStrategyChoice::All:
                return std::make_unique<AllMaskStrategy>();
            case settings::rigidbody::ParameterMaskStrategyChoice::Real:
                return std::make_unique<RealOnlyMaskStrategy>();
            case settings::rigidbody::ParameterMaskStrategyChoice::Symmetry:
                return std::make_unique<SymmetryOnlyMaskStrategy>();
            case settings::rigidbody::ParameterMaskStrategyChoice::Sequential:
                return std::make_unique<SequentialMaskStrategy>();
            case settings::rigidbody::ParameterMaskStrategyChoice::Random:
                return std::make_unique<RandomMaskStrategy>();
            case settings::rigidbody::ParameterMaskStrategyChoice::SequentialSymmetry:
                return std::make_unique<SequentialSymmetryMaskStrategy>();
            case settings::rigidbody::ParameterMaskStrategyChoice::SequentialReal:
                return std::make_unique<SequentialRealMaskStrategy>();
            default:
                throw except::unknown_argument("rigidbody::factory::create_mask_strategy: Unknown mask strategy. Did you forget to add it?");
        }
    }
}

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(observer_ptr<const Rigidbody> body) {
    return create_selection_strategy(body, settings::rigidbody::body_select_strategy, settings::rigidbody::parameter_mask_strategy);
}

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(observer_ptr<const Rigidbody> body, settings::rigidbody::BodySelectStrategyChoice choice) {
    return create_selection_strategy(body, choice, settings::rigidbody::parameter_mask_strategy);
}

std::unique_ptr<BodySelectStrategy> rigidbody::factory::create_selection_strategy(
    observer_ptr<const Rigidbody> body,
    settings::rigidbody::BodySelectStrategyChoice body_choice,
    settings::rigidbody::ParameterMaskStrategyChoice mask_choice)
{
    std::unique_ptr<BodySelectStrategy> strategy;
    switch (body_choice) {
        case settings::rigidbody::BodySelectStrategyChoice::RandomBodySelect:
            strategy = std::make_unique<RandomBodySelect>(body); break;
        case settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect:
            strategy = std::make_unique<RandomConstraintSelect>(body); break;
        case settings::rigidbody::BodySelectStrategyChoice::SequentialBodySelect:
            strategy = std::make_unique<SequentialBodySelect>(body); break;
        case settings::rigidbody::BodySelectStrategyChoice::SequentialConstraintSelect:
            strategy = std::make_unique<SequentialConstraintSelect>(body); break;
        default:
            throw except::unknown_argument("rigidbody::factory::create_selection_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
    strategy->set_mask_strategy(create_mask_strategy(mask_choice));
    return strategy;
}