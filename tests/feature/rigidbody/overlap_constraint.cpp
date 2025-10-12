#include "settings/MoleculeSettings.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    fixture() {
        settings::molecule::implicit_hydrogens = false;

        a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
        a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
        a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
        a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
        a5 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
        a6 = AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C);
        a7 = AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
        a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::NH);

        b1 = Body(std::vector{a1, a2});
        b2 = Body(std::vector{a3, a4});
        b3 = Body(std::vector{a5, a6});
        b4 = Body(std::vector{a7, a8});
        ap = {b1, b2, b3, b4};
    }

    AtomFF a1, a2, a3, a4, a5, a6, a7, a8;
    Body b1, b2, b3, b4;
    std::vector<Body> ap;
};

// We cannot directly test the constructor, only its effects on the evaluate() function.
// TEST_CASE_METHOD(fixture, "OverlapConstraint::OverlapConstraint") {
//     SECTION("Protein*") {
//         OverlapConstraint oc(&protein);
//         REQUIRE(oc.evaluate() == 0);
//     }
// }

TEST_CASE_METHOD(fixture, "OverlapConstraint::evaluate") {
    settings::general::verbose = false;

    Molecule protein(ap);
    SECTION("initializes to 0") {
        constraints::OverlapConstraint oc(&protein);
        REQUIRE(oc.evaluate() == 0);
    }

    SECTION("returns non-zero when moved closer") {
        constraints::OverlapConstraint oc(&protein);
        REQUIRE(oc.evaluate() == 0);
        
        protein.get_body(0).translate(Vector3<double>(2, 2, 1.5));
        REQUIRE(oc.evaluate() > 0);

        protein.get_body(0).translate(Vector3<double>(-2, -2, -1.5));
        REQUIRE(oc.evaluate() == 0);
    }
}