#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Constraint.h>
#include <rigidbody/ConstrainedFitter.h>
#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <rigidbody/transform/RigidTransform.h>
#include <data/Atom.h>

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
    Protein protein = Protein(ap);
};

TEST_CASE_METHOD(fixture, "constructor") {
    SECTION("Invalid constraints") {
        REQUIRE_THROWS(Constraint(&protein, a1, a2)); // same body
        REQUIRE_THROWS(Constraint(&protein, a6, a8)); // non-C
    }

    SECTION("Check construction") {
        Constraint c(&protein, a1, a3);

        REQUIRE(c.iatom1 == 0);
        REQUIRE(c.iatom2 == 0);
        REQUIRE(c.ibody1 == 0);
        REQUIRE(c.ibody2 == 1);
    }
}

TEST_CASE_METHOD(fixture, "connections") {
    RigidBody rigidbody(protein);
    
    RigidTransform transform(&rigidbody);
    std::shared_ptr<Constraint> constraint1 = std::make_shared<Constraint>(&protein, a1, a3);
    std::shared_ptr<Constraint> constraint2 = std::make_shared<Constraint>(&protein, a5, a7);
    rigidbody.add_constraint(constraint1);
    rigidbody.add_constraint(constraint2);

    std::vector<Body*> group = transform.get_connected(constraint1);
    REQUIRE(group.size() == 1);
    REQUIRE(*group[0] == b1);

    group = transform.get_connected(constraint2);
    REQUIRE(group.size() == 1);
    REQUIRE(*group[0] == b3);

    Constraint constraint3(&protein, a3, a5);
    rigidbody.add_constraint(constraint3);
    group = transform.get_connected(constraint3);
    REQUIRE(group.size() == 2);
    REQUIRE(*group[0] == b1);
    REQUIRE(*group[1] == b2);

    Constraint constraint4(&protein, a1, a7);
    rigidbody.add_constraint(constraint4);
    group = transform.get_connected(constraint3);
    REQUIRE(group.size() == 4);
}

TEST_CASE_METHOD(fixture, "affects_fitter") {
    setting::general::verbose = false;
    fitter::ConstrainedFitter fitter("test/files/2epe.dat", protein.get_histogram());
    double chi2 = fitter.fit()->fval;

    Constraint constraint(&protein, a1, a3);
    protein.body(0).translate(Vector3<double>(1, 0, 0));
    fitter.add_constraint(constraint);
    double chi2c = fitter.fit()->fval;

    CHECK(constraint.evaluate() > 0);
    REQUIRE_THAT(chi2c-chi2, Catch::Matchers::WithinAbs(constraint.evaluate(), 0.1));
}

TEST_CASE("simple_constraint_generation") {
    SECTION("simple") {
        int distance = setting::rigidbody::bond_distance;
        Atom a1 = Atom(Vector3<double>(0, 0, 0*distance), 1, "C", "C", 1);
        Atom a2 = Atom(Vector3<double>(0, 0, 1*distance), 1, "C", "C", 1);
        Atom a3 = Atom(Vector3<double>(0, 0, 2*distance), 1, "C", "C", 1);
        Atom a4 = Atom(Vector3<double>(0, 0, 3*distance), 1, "C", "C", 1);

        Body b1 = Body(std::vector<Atom>{a1});
        Body b2 = Body(std::vector<Atom>{a2});
        Body b3 = Body(std::vector<Atom>{a3});
        Body b4 = Body(std::vector<Atom>{a4});
        std::vector<Body> ap = {b1, b2, b3, b4};
        Protein protein = Protein(ap);

        RigidBody rigidbody(protein);
        rigidbody.generate_simple_constraints();
        REQUIRE(rigidbody.constraints.size() == 3);
    }

    SECTION("real data") {
        Protein protein = BodySplitter::split("data/LAR1-2/LAR1-2.pdb", {9, 99});
        RigidBody rigidbody(protein);
        rigidbody.generate_simple_constraints();

        REQUIRE(rigidbody.constraints.size() == 2);
        std::cout << rigidbody.constraints[0] << std::endl;
        std::cout << rigidbody.constraints[1] << std::endl;
    }
}