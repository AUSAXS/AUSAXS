#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <io/ExistingFile.h>
#include <settings/MoleculeSettings.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("VolumetricConstraints::generate") {
    settings::general::verbose = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric;
    settings::molecule::implicit_hydrogens = false;

    SECTION("simple") {
        double distance = settings::rigidbody::bond_distance;
        AtomFF a1({0,           0,           0}, form_factor::form_factor_t::C);
        AtomFF a2({0,           0,  1*distance}, form_factor::form_factor_t::C);
        AtomFF a3({0,           0, -1*distance}, form_factor::form_factor_t::C);
        AtomFF a4({0,  1*distance,           0}, form_factor::form_factor_t::C);
        AtomFF a5({0, -1*distance,           0}, form_factor::form_factor_t::C);
        AtomFF a6({ 1*distance, 0,           0}, form_factor::form_factor_t::C);
        AtomFF a7({-1*distance, 0,           0}, form_factor::form_factor_t::C);
        AtomFF a8({0,           0,  2*distance}, form_factor::form_factor_t::C);
        AtomFF a9({0,           0,  3*distance}, form_factor::form_factor_t::C);

        Body b1 = Body(std::vector<AtomFF>{a1});
        Body b2 = Body(std::vector<AtomFF>{a2});
        Body b3 = Body(std::vector<AtomFF>{a3});
        Body b4 = Body(std::vector<AtomFF>{a4});
        Body b5 = Body(std::vector<AtomFF>{a5});
        Body b6 = Body(std::vector<AtomFF>{a6});
        Body b7 = Body(std::vector<AtomFF>{a7});
        Body b8 = Body(std::vector<AtomFF>{a8});
        Body b9 = Body(std::vector<AtomFF>{a9});
        std::vector<Body> ap = {b1, b2, b3, b4, b5, b6, b7, b8, b9};
        rigidbody::Rigidbody rigidbody(Molecule{ap});
        REQUIRE(rigidbody.constraints->discoverable_constraints.size() == 8);
    }
}