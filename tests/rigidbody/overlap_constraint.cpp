#include "settings/MoleculeSettings.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace data;
using namespace data::record;
using namespace rigidbody;

struct fixture {
    fixture() {
        settings::molecule::use_effective_charge = false;
        settings::molecule::implicit_hydrogens = false;

        a1 = Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1);
        a2 = Atom(Vector3<double>(-1,  1, -1), 1, constants::atom_t::C, "C", 1);
        a3 = Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1);
        a4 = Atom(Vector3<double>(-1,  1,  1), 1, constants::atom_t::C, "C", 1);
        a5 = Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1);
        a6 = Atom(Vector3<double>( 1,  1, -1), 1, constants::atom_t::C, "C", 1);
        a7 = Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1);
        a8 = Atom(Vector3<double>( 1,  1,  1), 1, constants::atom_t::He, "He", 1);

        b1 = Body(std::vector<Atom>{a1, a2});
        b2 = Body(std::vector<Atom>{a3, a4});
        b3 = Body(std::vector<Atom>{a5, a6});
        b4 = Body(std::vector<Atom>{a7, a8});
        ap = {b1, b2, b3, b4};
    }

    Atom a1, a2, a3, a4, a5, a6, a7, a8;
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