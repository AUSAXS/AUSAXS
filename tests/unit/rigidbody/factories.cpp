#include <catch2/catch_test_macros.hpp>

#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/parameters/decay/LinearDecay.h>
#include <rigidbody/parameters/decay/ExponentialDecay.h>
#include <rigidbody/parameters/decay/NoDecay.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("DecayFactory::create_decay_strategy") {
    SECTION("Linear") {
        auto decay = factory::create_decay_strategy(100, settings::rigidbody::DecayStrategyChoice::Linear);
        REQUIRE(dynamic_cast<parameter::decay::LinearDecay*>(decay.get()) != nullptr);
    }

    SECTION("Exponential") {
        auto decay = factory::create_decay_strategy(100, settings::rigidbody::DecayStrategyChoice::Exponential);
        REQUIRE(dynamic_cast<parameter::decay::ExponentialDecay*>(decay.get()) != nullptr);
    }

    SECTION("None") {
        auto decay = factory::create_decay_strategy(100, settings::rigidbody::DecayStrategyChoice::None);
        REQUIRE(dynamic_cast<parameter::decay::NoDecay*>(decay.get()) != nullptr);
    }
}

TEST_CASE("BodySelectFactory::create_selection_strategy") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{a1}), Body(std::vector{a2})}});

    SECTION("RandomBodySelect") {
        auto strat = factory::create_selection_strategy(&rb, settings::rigidbody::BodySelectStrategyChoice::RandomBodySelect);
        REQUIRE(dynamic_cast<selection::RandomBodySelect*>(strat.get()) != nullptr);
    }

    SECTION("SequentialBodySelect") {
        auto strat = factory::create_selection_strategy(&rb, settings::rigidbody::BodySelectStrategyChoice::SequentialBodySelect);
        REQUIRE(dynamic_cast<selection::SequentialBodySelect*>(strat.get()) != nullptr);
    }
}

TEST_CASE("TransformFactory::create_transform_strategy") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{a1}), Body(std::vector{a2})}});

    SECTION("SingleTransform") {
        auto strat = factory::create_transform_strategy(&rb, settings::rigidbody::TransformationStrategyChoice::SingleTransform);
        REQUIRE(dynamic_cast<transform::SingleTransform*>(strat.get()) != nullptr);
    }

    SECTION("RigidTransform") {
        auto strat = factory::create_transform_strategy(&rb, settings::rigidbody::TransformationStrategyChoice::RigidTransform);
        REQUIRE(dynamic_cast<transform::RigidTransform*>(strat.get()) != nullptr);
    }
}

TEST_CASE("ControllerFactory::create_controller") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{a1}), Body(std::vector{a2})}});

    SECTION("Classic") {
        auto ctrl = factory::create_controller(&rb, settings::rigidbody::ControllerChoice::Classic);
        REQUIRE(dynamic_cast<controller::SimpleController*>(ctrl.get()) != nullptr);
    }
}
