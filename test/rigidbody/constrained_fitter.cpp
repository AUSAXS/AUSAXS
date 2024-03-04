#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/constraints/ConstrainedFitter.h>
#include <hist//intensity_calculator/ICompositeDistanceHistogram.h>
#include <fitter/HydrationFitter.h>
#include <rigidbody/RigidBody.h>
#include <math/Vector3.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <settings/All.h>

using namespace data;
using namespace data::record;
using namespace rigidbody;

struct fixture {
    Atom a1 = Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1);
    Atom a2 = Atom(Vector3<double>(-1,  1, -1), 1, constants::atom_t::C, "C", 1);
    Atom a3 = Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1);
    Atom a4 = Atom(Vector3<double>(-1,  1,  1), 1, constants::atom_t::C, "C", 1);
    Atom a5 = Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1);
    Atom a6 = Atom(Vector3<double>( 1,  1, -1), 1, constants::atom_t::C, "C", 1);
    Atom a7 = Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1);
    Atom a8 = Atom(Vector3<double>( 1,  1,  1), 1, constants::atom_t::He, "He", 1);

    Body b1 = Body(std::vector<Atom>{a1, a2});
    Body b2 = Body(std::vector<Atom>{a3, a4});
    Body b3 = Body(std::vector<Atom>{a5, a6});
    Body b4 = Body(std::vector<Atom>{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

TEST_CASE_METHOD(fixture, "ConstrainedFitter::constraint_manager") {
    settings::general::verbose = false;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    RigidBody protein(ap);

    fitter::ConstrainedFitter<fitter::HydrationFitter> fitter("test/files/2epe.dat", protein.get_histogram());
    CHECK(fitter.get_constraint_manager() == nullptr);
    fitter.set_constraint_manager(protein.get_constraint_manager());
    CHECK(fitter.get_constraint_manager() == protein.get_constraint_manager().get());
}

TEST_CASE_METHOD(fixture, "ConstrainedFitter::chi2") {
    settings::general::verbose = false;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    RigidBody protein(ap);

    fitter::ConstrainedFitter<fitter::HydrationFitter> fitter("test/files/2epe.dat", protein.get_histogram());
    fitter.set_constraint_manager(protein.get_constraint_manager());
    double chi2 = fitter.fit()->fval;

    constraints::DistanceConstraint constraint(&protein, a1, a3);
    protein.get_body(0).translate(Vector3<double>(1, 0, 0));
    fitter.get_constraint_manager()->add_constraint(std::move(constraint));
    double chi2c = fitter.fit()->fval;

    CHECK(constraint.evaluate() > 0);
    REQUIRE_THAT(chi2c-chi2, Catch::Matchers::WithinAbs(constraint.evaluate(), 0.1));
}