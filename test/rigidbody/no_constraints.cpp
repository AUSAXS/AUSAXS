#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Atom.h>
#include <data/Body.h>
#include <data/BodySplitter.h>
#include <rigidbody/RigidBody.h>
#include <io/ExistingFile.h>
#include <settings/ProteinSettings.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>

TEST_CASE("LinearConstraints::generate") {
    settings::general::verbose = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::protein::use_effective_charge = false;

    SECTION("simple") {
        int distance = settings::rigidbody::bond_distance;
        Atom a1 = Atom(Vector3<double>(0, 0, 0*distance), 1, "C", "C", 1);
        Atom a2 = Atom(Vector3<double>(0, 0, 1*distance), 1, "C", "C", 1);
        Atom a3 = Atom(Vector3<double>(0, 0, 2*distance), 1, "C", "C", 1);
        Atom a4 = Atom(Vector3<double>(0, 0, 3*distance), 1, "C", "C", 1);

        Body b1 = Body(std::vector<Atom>{a1});
        Body b2 = Body(std::vector<Atom>{a2});
        Body b3 = Body(std::vector<Atom>{a3});
        Body b4 = Body(std::vector<Atom>{a4});
        std::vector<Body> ap = {b1, b2, b3, b4};
        rigidbody::RigidBody rigidbody(ap);
        REQUIRE(rigidbody.get_constraint_manager()->distance_constraints.size() == 0);
    }

    SECTION("real data") {
        rigidbody::RigidBody rigidbody = rigidbody::BodySplitter::split("test/files/LAR1-2.pdb", {9, 99});
        REQUIRE(rigidbody.get_constraint_manager()->distance_constraints.size() == 0);
    }
}