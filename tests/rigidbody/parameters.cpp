#include <catch2/catch_test_macros.hpp>

#include <rigidbody/parameters/Parameters.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/Body.h>
#include <settings/GeneralSettings.h>
#include <settings/MoleculeSettings.h>

using namespace ausaxs;
using namespace data;
using namespace data::record;
using namespace rigidbody;

struct fixture {
    Atom a1 = Atom(1, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1, -1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a2 = Atom(2, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1,  1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a3 = Atom(3, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1, -1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a4 = Atom(4, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1,  1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a5 = Atom(5, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1, -1,  1), 1, 0, constants::atom_t::C, "0");
    Atom a6 = Atom(6, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1,  1,  1), 1, 0, constants::atom_t::C, "0");
    Atom a7 = Atom(7, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1, -1,  1), 1, 0, constants::atom_t::C, "0");
    Atom a8 = Atom(8, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1,  1,  1), 1, 0, constants::atom_t::C, "0");
    
    std::vector<Atom> b1 = {a1, a2};
    std::vector<Atom> b2 = {a3, a4};
    std::vector<Atom> b3 = {a5, a6};
    std::vector<Atom> b4 = {a7, a8};
    std::vector<Body> bodies = {Body(b1), Body(b2), Body(b3), Body(b4)};
};

TEST_CASE("Parameters::Parameter") {
    SECTION("default") {
        rigidbody::parameter::Parameter p;
        CHECK(p.dr == Vector3<double>(0, 0, 0));
        CHECK(p.alpha == 0);
        CHECK(p.beta == 0);
        CHECK(p.gamma == 0);
    }

    SECTION("Vector3<double>&, double, double, double") {
        Vector3<double> dx(1, 2, 3);
        double alpha = 4, beta = 5, gamma = 6;
        rigidbody::parameter::Parameter p(dx, alpha, beta, gamma);
        CHECK(p.dr == dx);
        CHECK(p.alpha == alpha);
        CHECK(p.beta == beta);
        CHECK(p.gamma == gamma);
    }
}

TEST_CASE_METHOD(fixture, "Parameters::Parameters") {
    settings::general::verbose = false;

    SECTION("Protein*") {
        Molecule protein(bodies);
        parameter::Parameters params(&protein);
        CHECK(params.params.size() == 4);
        CHECK(params.id_to_index.size() == 4);
    }
}    

TEST_CASE_METHOD(fixture, "Parameters::update") {
    // also tests Parameters::get
    settings::general::verbose = false;

    SECTION("unsigned int, const Parameter&") {
        Molecule protein(bodies);
        parameter::Parameters params(&protein);
        parameter::Parameter p(Vector3<double>(1, 2, 3), 4, 5, 6);

        auto id = protein.get_body(1).get_id();
        params.update(id, p);
        CHECK(params.get(id).dr == p.dr);
        CHECK(params.get(id).alpha == p.alpha);
        CHECK(params.get(id).beta == p.beta);
        CHECK(params.get(id).gamma == p.gamma);
    }

    SECTION("unsigned int, Vector3<double>, double, double, double") {
        Molecule protein(bodies);
        parameter::Parameters params(&protein);
        Vector3<double> dx(1, 2, 3);
        double alpha = 4, beta = 5, gamma = 6;
        
        auto id = protein.get_body(1).get_id();
        params.update(id, dx, alpha, beta, gamma);
        CHECK(params.get(id).dr == dx);
        CHECK(params.get(id).alpha == alpha);
        CHECK(params.get(id).beta == beta);
        CHECK(params.get(id).gamma == gamma);
    }
}