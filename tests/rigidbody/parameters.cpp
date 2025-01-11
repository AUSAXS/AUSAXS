#include <catch2/catch_test_macros.hpp>

#include <rigidbody/parameters/Parameters.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/GeneralSettings.h>
#include <settings/MoleculeSettings.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a5 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a6 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    AtomFF a7 = AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::C);
    
    std::vector<AtomFF> b1 = {a1, a2};
    std::vector<AtomFF> b2 = {a3, a4};
    std::vector<AtomFF> b3 = {a5, a6};
    std::vector<AtomFF> b4 = {a7, a8};
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

        auto id = protein.get_body(1).get_uid();
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
        
        auto id = protein.get_body(1).get_uid();
        params.update(id, dx, alpha, beta, gamma);
        CHECK(params.get(id).dr == dx);
        CHECK(params.get(id).alpha == alpha);
        CHECK(params.get(id).beta == beta);
        CHECK(params.get(id).gamma == gamma);
    }
}