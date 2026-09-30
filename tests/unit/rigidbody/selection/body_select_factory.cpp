#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("BodySelectFactory::create_selection_strategy") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

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
