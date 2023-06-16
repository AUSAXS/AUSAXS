#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/Protein.h>
#include <data/Atom.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace rigidbody;

struct fixture {
    Atom a1 = Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1);
    Atom a2 = Atom(Vector3<double>(-1,  1, -1), 1, "C", "C", 1);
    Atom a3 = Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1);
    Atom a4 = Atom(Vector3<double>(-1,  1,  1), 1, "C", "C", 1);
    Atom a5 = Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1);
    Atom a6 = Atom(Vector3<double>( 1,  1, -1), 1, "C", "C", 1);
    Atom a7 = Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1);
    Atom a8 = Atom(Vector3<double>( 1,  1,  1), 1, "He", "He", 1);

    Body b1 = Body(std::vector<Atom>{a1, a2});
    Body b2 = Body(std::vector<Atom>{a3, a4});
    Body b3 = Body(std::vector<Atom>{a5, a6});
    Body b4 = Body(std::vector<Atom>{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

// We cannot directly test the constructor, only its effects on the evaluate() function.
// TEST_CASE_METHOD(fixture, "OverlapConstraint::OverlapConstraint") {
//     SECTION("Protein*") {
//         OverlapConstraint oc(&protein);
//         REQUIRE(oc.evaluate() == 0);
//     }
// }

TEST_CASE_METHOD(fixture, "OverlapConstraint::evaluate") {
    settings::protein::use_effective_charge = false;
    settings::axes::distance_bin_width = 0.1;
    Protein protein(ap);
    SECTION("initializes to 0") {
        OverlapConstraint oc(&protein);
        REQUIRE(oc.evaluate() == 0);
    }

    SECTION("returns non-zero when moved closer") {
        OverlapConstraint oc(&protein);
        REQUIRE(oc.evaluate() == 0);
        
        protein.get_body(0).translate(Vector3<double>(2, 2, 1.5));
        REQUIRE(oc.evaluate() > 0);

        protein.get_body(0).translate(Vector3<double>(-2, -2, -1.5));
        REQUIRE(oc.evaluate() == 0);
    }
    settings::axes::distance_bin_width = 1;
}