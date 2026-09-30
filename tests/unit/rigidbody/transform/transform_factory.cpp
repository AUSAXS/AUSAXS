#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformFactory.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("TransformFactory::create_transform_strategy") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{a1}), Body(std::vector{a2})}});

    SECTION("SingleTransform") {
        auto strat = factory::create_transform_strategy(&rb, settings::rigidbody::TransformationStrategyChoice::SingleTransform);
        REQUIRE(dynamic_cast<rigidbody::transform::SingleTransform*>(strat.get()) != nullptr);
    }

    SECTION("RigidTransform") {
        auto strat = factory::create_transform_strategy(&rb, settings::rigidbody::TransformationStrategyChoice::RigidTransform);
        REQUIRE(dynamic_cast<rigidbody::transform::RigidTransform*>(strat.get()) != nullptr);
    }
}
