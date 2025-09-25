#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <hist//intensity_calculator/ICompositeDistanceHistogram.h>
#include <fitter/SmartFitter.h>
#include <rigidbody/Rigidbody.h>
#include <math/Vector3.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    AtomFF a5 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a6 = AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a7 = AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::NH);

    Body b1 = Body(std::vector<AtomFF>{a1, a2});
    Body b2 = Body(std::vector<AtomFF>{a3, a4});
    Body b3 = Body(std::vector<AtomFF>{a5, a6});
    Body b4 = Body(std::vector<AtomFF>{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

TEST_CASE_METHOD(fixture, "ConstrainedFitter::constraint_manager") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    Rigidbody protein(Molecule{ap});

    fitter::ConstrainedFitter fitter(protein.constraints.get(), SimpleDataset("tests/files/2epe.dat"), protein.molecule.get_histogram());
    CHECK(fitter.constraints == protein.constraints.get());
}

TEST_CASE_METHOD(fixture, "ConstrainedFitter::chi2") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None; // make sure there's no other distance constraints
    Rigidbody protein(Molecule{ap});

    fitter::ConstrainedFitter fitter(protein.constraints.get(), SimpleDataset("tests/files/2epe.dat"), protein.molecule.get_histogram());
    double chi2 = fitter.fit()->fval;

    constraints::DistanceConstraint constraint(&protein.molecule, a1, a3);
    protein.molecule.get_body(0).translate(Vector3<double>(-1, 0, 0));
    fitter.constraints->add_constraint(std::move(constraint));
    double chi2c = fitter.fit()->fval;

    CHECK(constraint.evaluate() > 0);
    REQUIRE_THAT(chi2c-chi2, Catch::Matchers::WithinAbs(constraint.evaluate(), 0.1));
}