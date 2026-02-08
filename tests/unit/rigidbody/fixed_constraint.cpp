#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/FixedConstraint.h>
#include <data/Molecule.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;

TEST_CASE("FixedConstraint::construction") {
    data::Molecule molecule("tests/files/SASDJG5.pdb");
    
    SECTION("construct with two different bodies") {
        FixedConstraint constraint(&molecule, 0, 1);
        
        CHECK(constraint.ibody1 == 0);
        CHECK(constraint.ibody2 == 1);
        CHECK(constraint.protein == &molecule);
    }

    SECTION("throws on same body constraint") {
        CHECK_THROWS(FixedConstraint(&molecule, 0, 0));
    }
}

TEST_CASE("FixedConstraint::evaluate") {
    data::Molecule molecule("tests/files/SASDJG5.pdb");
    
    SECTION("always returns zero") {
        FixedConstraint constraint(&molecule, 0, 1);
        CHECK(constraint.evaluate() == 0.0);
    }
}

TEST_CASE("FixedConstraint::get_body") {
    data::Molecule molecule("tests/files/SASDJG5.pdb");
    FixedConstraint constraint(&molecule, 0, 1);
    
    SECTION("get_body1 returns first body") {
        const auto& body1 = constraint.get_body1();
        CHECK(&body1 == &molecule.get_body(0));
    }

    SECTION("get_body2 returns second body") {
        const auto& body2 = constraint.get_body2();
        CHECK(&body2 == &molecule.get_body(1));
    }

    SECTION("non-const get_body1 returns first body") {
        FixedConstraint non_const_constraint(&molecule, 0, 1);
        auto& body1 = non_const_constraint.get_body1();
        CHECK(&body1 == &molecule.get_body(0));
    }

    SECTION("non-const get_body2 returns second body") {
        FixedConstraint non_const_constraint(&molecule, 0, 1);
        auto& body2 = non_const_constraint.get_body2();
        CHECK(&body2 == &molecule.get_body(1));
    }
}
