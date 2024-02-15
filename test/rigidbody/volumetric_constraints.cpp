#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/record/Atom.h>
#include <data/Body.h>
#include <data/BodySplitter.h>
#include <rigidbody/RigidBody.h>
#include <io/ExistingFile.h>
#include <settings/MoleculeSettings.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>

using namespace data;
using namespace data::record;

TEST_CASE("VolumetricConstraints::generate") {
    settings::general::verbose = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric;
    settings::molecule::use_effective_charge = false;

    SECTION("simple") {
        int distance = settings::rigidbody::bond_distance;
        Atom a1 = Atom(Vector3<double>(0,  0,  0), 1, constants::atom_t::C, "C", 1);
        Atom a2 = Atom(Vector3<double>(0,  0,  1*distance), 1, constants::atom_t::C, "C", 1);
        Atom a3 = Atom(Vector3<double>(0,  0, -1*distance), 1, constants::atom_t::C, "C", 1);
        Atom a4 = Atom(Vector3<double>(0,  1*distance,  0), 1, constants::atom_t::C, "C", 1);
        Atom a5 = Atom(Vector3<double>(0, -1*distance,  0), 1, constants::atom_t::C, "C", 1);
        Atom a6 = Atom(Vector3<double>( 1*distance, 0,  0), 1, constants::atom_t::C, "C", 1);
        Atom a7 = Atom(Vector3<double>(-1*distance, 0,  0), 1, constants::atom_t::C, "C", 1);
        Atom a8 = Atom(Vector3<double>(0,  0,  2*distance), 1, constants::atom_t::C, "C", 1);
        Atom a9 = Atom(Vector3<double>(0,  0,  3*distance), 1, constants::atom_t::C, "C", 1);

        Body b1 = Body(std::vector<Atom>{a1});
        Body b2 = Body(std::vector<Atom>{a2});
        Body b3 = Body(std::vector<Atom>{a3});
        Body b4 = Body(std::vector<Atom>{a4});
        Body b5 = Body(std::vector<Atom>{a5});
        Body b6 = Body(std::vector<Atom>{a6});
        Body b7 = Body(std::vector<Atom>{a7});
        Body b8 = Body(std::vector<Atom>{a8});
        Body b9 = Body(std::vector<Atom>{a9});
        std::vector<Body> ap = {b1, b2, b3, b4, b5, b6, b7, b8, b9};
        rigidbody::RigidBody rigidbody(ap);
        REQUIRE(rigidbody.get_constraint_manager()->distance_constraints.size() == 8);
    }
}