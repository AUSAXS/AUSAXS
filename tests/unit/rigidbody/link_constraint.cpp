#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/LinkConstraint.h>
#include <data/Molecule.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;

TEST_CASE("LinkConstraint::construction") {
    data::Molecule molecule("tests/files/SASDJG5.pdb");
    
    SECTION("construct with two different bodies") {
        LinkConstraint constraint(&molecule, 0, 1);
        
        CHECK(constraint.ibody1 == 0);
        CHECK(constraint.ibody2 == 1);
        CHECK(constraint.protein == &molecule);
    }

    SECTION("throws on same body constraint") {
        CHECK_THROWS(LinkConstraint(&molecule, 0, 0));
    }
}

TEST_CASE("LinkConstraint::evaluate") {
    data::Molecule molecule("tests/files/SASDJG5.pdb");
    
    SECTION("always returns zero") {
        LinkConstraint constraint(&molecule, 0, 1);
        CHECK(constraint.evaluate() == 0.0);
    }
}

TEST_CASE("LinkConstraint::get_body") {
    data::Molecule molecule("tests/files/SASDJG5.pdb");
    LinkConstraint constraint(&molecule, 0, 1);
    
    SECTION("get_body1 returns first body") {
        const auto& body1 = constraint.get_body1();
        CHECK(&body1 == &molecule.get_body(0));
    }

    SECTION("get_body2 returns second body") {
        const auto& body2 = constraint.get_body2();
        CHECK(&body2 == &molecule.get_body(1));
    }

    SECTION("non-const get_body1 returns first body") {
        LinkConstraint non_const_constraint(&molecule, 0, 1);
        auto& body1 = non_const_constraint.get_body1();
        CHECK(&body1 == &molecule.get_body(0));
    }

    SECTION("non-const get_body2 returns second body") {
        LinkConstraint non_const_constraint(&molecule, 0, 1);
        auto& body2 = non_const_constraint.get_body2();
        CHECK(&body2 == &molecule.get_body(1));
    }
}
