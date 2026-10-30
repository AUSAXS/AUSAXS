#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/GroupFormFactor.h>
#include <form_factor/neutron/NeutronFormFactor.h>

#include <cmath>

using namespace ausaxs;
using namespace form_factor;
using Catch::Matchers::WithinAbs;

TEST_CASE("GroupFormFactor::evaluate") {
    // CH2 with neutron scattering lengths, which do not depend on q
    double bc = 6.646, bh = -3.739, dxh = 1.083, dhh = 1.768;
    GroupFormFactor<neutron::AtomicFormFactor> ff(bc, bh, {.hydrogens = 2, .d_xh = dxh, .d_hh = dhh});
    auto sinc = [] (double x) {return std::sin(x)/x;};

    SECTION("at q = 0") {
        CHECK_THAT(ff.evaluate(0), WithinAbs(bc + 2*bh, 1e-12));
        CHECK_THAT(ff.I0(), WithinAbs(bc + 2*bh, 1e-12));
    }

    SECTION("at non-zero q") {
        CHECK_THAT(ff.evaluate(1), WithinAbs(bc + 2*bh*sinc(dxh), 1e-12));
    }

    SECTION("self-term equals the squared amplitude at q = 0") {
        CHECK_THAT(ff.evaluate_self(0), WithinAbs(ff.I0()*ff.I0(), 1e-12));
    }

    SECTION("self-term at non-zero q") {
        double expected = bc*bc + 2*bh*bh + 4*bc*bh*sinc(dxh) + 2*bh*bh*sinc(dhh);
        CHECK_THAT(ff.evaluate_self(1), WithinAbs(expected, 1e-12));
        CHECK(ff.evaluate_self(1) > 10*std::pow(ff.evaluate(1), 2));
    }

    SECTION("a lone atom is its own form factor") {
        GroupFormFactor<neutron::AtomicFormFactor> atom(bc, bh, group::atom);
        CHECK_THAT(atom.evaluate(1), WithinAbs(bc, 1e-12));
        CHECK_THAT(atom.evaluate_self(1), WithinAbs(bc*bc, 1e-12));
        CHECK_THAT(atom.I0(), WithinAbs(bc, 1e-12));
    }

    SECTION("compile-time evaluation matches run-time evaluation") {
        constexpr GroupFormFactor<neutron::AtomicFormFactor> cff(6.646, -3.739, {.hydrogens = 2, .d_xh = 1.083, .d_hh = 1.768});
        constexpr double val = cff.evaluate(0.5);
        constexpr double val_self = cff.evaluate_self(0.5);
        CHECK_THAT(val, WithinAbs(ff.evaluate(0.5), 1e-12));
        CHECK_THAT(val_self, WithinAbs(ff.evaluate_self(0.5), 1e-12));
    }
}

TEST_CASE("group::hh_distance") {
    // two hydrogens at right angles are sqrt(2) bond lengths apart
    CHECK_THAT(group::hh_distance(1, 90), WithinAbs(std::sqrt(2), 1e-12));
    CHECK_THAT(group::hh_distance(group::bond::CH, group::angle::tetrahedral), WithinAbs(1.768, 1e-3));
}
