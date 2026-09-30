#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <settings/All.h>

#include <support/rb_metadata.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody::constraints;

struct fixture {
    fixture() {
        settings::molecule::implicit_hydrogens = false;
    }

    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    AtomFF a5 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a6 = AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a7 = AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::NH);

    Body b1 = Body(std::vector{a1, a2});
    Body b2 = Body(std::vector{a3, a4});
    Body b3 = Body(std::vector{a5, a6});
    Body b4 = Body(std::vector{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

TEST_CASE_METHOD(fixture, "OverlapConstraint::evaluate") {
    Molecule mol(ap);
    test::mark_backbone_carbons(mol);

    OverlapConstraint o(&mol);
    // initialize() runs in ctor; evaluate should be non-negative
    double v0 = o.evaluate();
    CHECK(v0 >= 0);
    mol.get_body(0).translate(Vector3<double>(0,0,1));
    double v1 = o.evaluate();
    CHECK(v1 >= 0);
    // values may change after translation
    CHECK(v1 != v0);
}
