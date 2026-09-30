#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/controller/SimpleController.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("ControllerFactory::create_controller") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{a1}), Body(std::vector{a2})}});

    SECTION("Classic") {
        auto ctrl = factory::create_controller(&rb, settings::rigidbody::ControllerChoice::Classic);
        REQUIRE(dynamic_cast<controller::SimpleController*>(ctrl.get()) != nullptr);
    }
}
