#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/neutron/NeutronFormFactor.h>

#include <cmath>

using namespace ausaxs;
using namespace form_factor;
using Catch::Matchers::WithinAbs;

TEST_CASE("NeutronFormFactor::evaluate") {
    // A heavy atom with two hydrogens at some arbitrary distance. 
    constexpr double bc = 1, dxh = 3, dhh = 5;
    constexpr double bh = neutron::scattering_length::H;
    constexpr neutron::FormFactor ff({.b_heavy = bc, .hydrogens = 2, .d_xh = dxh, .d_hh = dhh});
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
        CHECK(std::abs(ff.evaluate_self(1) - std::pow(ff.evaluate(1), 2)) > 1); // a group is not a point scatterer
    }

    SECTION("compile-time evaluation matches run-time evaluation") {
        constexpr double val = ff.evaluate(0.5);
        constexpr double val_self = ff.evaluate_self(0.5);
        CHECK_THAT(val, WithinAbs(ff.evaluate(0.5), 1e-12));
        CHECK_THAT(val_self, WithinAbs(ff.evaluate_self(0.5), 1e-12));
    }
}

TEST_CASE("NeutronFormFactor::get") {
    using namespace neutron;

    SECTION("types without hydrogens are q-independent") {
        for (auto type : {form_factor_t::C, form_factor_t::N, form_factor_t::O, form_factor_t::S, form_factor_t::OTHER}) {
            const auto& ff = get(type);
            CHECK_THAT(ff.evaluate(1), WithinAbs(ff.I0(), 1e-12));
        }
    }

    SECTION("types with hydrogens vary with q") {
        for (auto type : {form_factor_t::CH, form_factor_t::CH2, form_factor_t::CH3, form_factor_t::NH, form_factor_t::NH2, form_factor_t::NH3, form_factor_t::OH, form_factor_t::SH}) {
            const auto& ff = get(type);
            CHECK(std::abs(ff.evaluate(1) - ff.I0()) > 0.1);
        }
    }

    SECTION("self-term equals the squared amplitude for single atoms") {
        for (auto type : {form_factor_t::H, form_factor_t::C, form_factor_t::N, form_factor_t::O, form_factor_t::S, form_factor_t::OTHER}) {
            const auto& ff = get(type);
            CHECK_THAT(ff.evaluate_self(1), WithinAbs(ff.evaluate(1)*ff.evaluate(1), 1e-12));
        }
    }

    SECTION("self-term is non-negative") {
        for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
            for (double q : {0., 0.25, 0.5, 1., 2., 5.}) {
                CHECK(get(static_cast<form_factor_t>(i)).evaluate_self(q) >= 0);
            }
        }
    }

    SECTION("excluded volume is not defined") {
        CHECK_THROWS(get(form_factor_t::EXCLUDED_VOLUME));
    }
}

TEST_CASE("neutron::solvent_density") {
    // the tabulated SLDs are -0.56e-6 and 6.36e-6 Å^-2, and 1 fm/Å^3 = 1e-5 Å^-2
    CHECK_THAT(neutron::solvent_density(0), WithinAbs(-0.056, 0.001));
    CHECK_THAT(neutron::solvent_density(1), WithinAbs(0.636, 0.001));
    CHECK_THAT(neutron::solvent_density(0.25), Catch::Matchers::WithinRel(0.75*neutron::solvent_density(0) + 0.25*neutron::solvent_density(1), 1e-12));
    CHECK_THROWS(neutron::solvent_density(-0.1));
    CHECK_THROWS(neutron::solvent_density(1.1));
}
