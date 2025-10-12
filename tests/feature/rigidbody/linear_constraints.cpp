#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Body.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/RigidBody.h>
#include <io/ExistingFile.h>
#include <settings/RigidBodySettings.h>
#include <settings/MoleculeSettings.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("LinearConstraints::generate") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    SECTION("simple") {
        double distance = settings::rigidbody::bond_distance;
        AtomFF a1({0, 0, 0*distance}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 1*distance}, form_factor::form_factor_t::C);
        AtomFF a3({0, 0, 2*distance}, form_factor::form_factor_t::C);
        AtomFF a4({0, 0, 3*distance}, form_factor::form_factor_t::C);

        Body b1 = Body(std::vector<AtomFF>{a1});
        Body b2 = Body(std::vector<AtomFF>{a2});
        Body b3 = Body(std::vector<AtomFF>{a3});
        Body b4 = Body(std::vector<AtomFF>{a4});
        std::vector<Body> ap = {b1, b2, b3, b4};
        rigidbody::RigidBody rigidbody(ap);
        REQUIRE(rigidbody.get_constraint_manager()->distance_constraints.size() == 3);
    }

    SECTION("real data") {
        rigidbody::RigidBody rigidbody = rigidbody::BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        REQUIRE(rigidbody.get_constraint_manager()->distance_constraints.size() == 2);
    }
}