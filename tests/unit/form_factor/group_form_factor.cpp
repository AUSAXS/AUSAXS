#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/GroupFormFactor.h>
#include <form_factor/neutron/NeutronFormFactor.h>

#include <cmath>

using namespace ausaxs;
using namespace form_factor;
using Catch::Matchers::WithinAbs;

TEST_CASE("GroupFormFactor::evaluate") {
    // arbitrary vals
    constexpr double bc = 3, bh = -2;
    constexpr group::Geometry geometry {.hydrogens = 4, .d_xh = 1.5, .d_hh = 2.5};
    const auto& [n, dxh, dhh] = geometry;
    GroupFormFactor<neutron::AtomicFormFactor> ff(bc, bh, geometry);
    auto sinc = [] (double x) {return std::sin(x)/x;};

    SECTION("at q = 0") {
        CHECK_THAT(ff.evaluate(0), WithinAbs(bc + n*bh, 1e-12));
        CHECK_THAT(ff.I0(), WithinAbs(bc + n*bh, 1e-12));
    }

    SECTION("at non-zero q") {
        CHECK_THAT(ff.evaluate(1), WithinAbs(bc + n*bh*sinc(dxh), 1e-12));
    }

    SECTION("self-term equals the squared amplitude at q = 0") {
        CHECK_THAT(ff.evaluate_self(0), WithinAbs(ff.I0()*ff.I0(), 1e-12));
    }

    SECTION("self-term at non-zero q") {
        double expected = bc*bc + n*bh*bh + 2*n*bc*bh*sinc(dxh) + n*(n-1)*bh*bh*sinc(dhh);
        CHECK_THAT(ff.evaluate_self(1), WithinAbs(expected, 1e-12));
    }

    SECTION("a lone atom is its own form factor") {
        GroupFormFactor<neutron::AtomicFormFactor> atom(bc, bh, group::atom);
        CHECK_THAT(atom.evaluate(1), WithinAbs(bc, 1e-12));
        CHECK_THAT(atom.evaluate_self(1), WithinAbs(bc*bc, 1e-12));
        CHECK_THAT(atom.I0(), WithinAbs(bc, 1e-12));
    }

    SECTION("compile-time evaluation matches run-time evaluation") {
        constexpr GroupFormFactor<neutron::AtomicFormFactor> cff(bc, bh, geometry);
        constexpr double val = cff.evaluate(0.5);
        constexpr double val_self = cff.evaluate_self(0.5);
        CHECK_THAT(val, WithinAbs(ff.evaluate(0.5), 1e-12));
        CHECK_THAT(val_self, WithinAbs(ff.evaluate_self(0.5), 1e-12));
    }
}

TEST_CASE("group::hh_distance") {
    // two hydrogens at right angles are sqrt(2) bond lengths apart
    CHECK_THAT(group::hh_distance(1, 90), WithinAbs(std::sqrt(2), 1e-12));

    // sin(tetrahedral/2) = sqrt(2/3), so the H-H distance is sqrt(8/3) bond lengths
    CHECK_THAT(group::hh_distance(2, group::angle::tetrahedral), WithinAbs(2*std::sqrt(8.0/3), 1e-6));
}
