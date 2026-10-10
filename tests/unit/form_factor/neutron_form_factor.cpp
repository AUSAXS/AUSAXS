#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/NeutronFormFactor.h>

#include <cmath>

using namespace ausaxs;
using namespace form_factor;
using Catch::Matchers::WithinAbs;

TEST_CASE("NeutronFormFactor::evaluate") {
    // CH2 with protium
    double bc = 6.646, bh = -3.739, dxh = 1.083, dhh = 1.768;
    neutron::FormFactor ff(bc, bh, 2, dxh, dhh);
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

    SECTION("compile-time evaluation matches run-time evaluation") {
        constexpr neutron::FormFactor cff(6.646, -3.739, 2, 1.083, 1.768);
        constexpr double val = cff.evaluate(0.5);
        constexpr double val_self = cff.evaluate_self(0.5);
        CHECK_THAT(val, WithinAbs(ff.evaluate(0.5), 1e-12));
        CHECK_THAT(val_self, WithinAbs(ff.evaluate_self(0.5), 1e-12));
    }
}

TEST_CASE("NeutronFormFactor::lookup") {
    using namespace neutron;

    SECTION("q = 0 values are the summed scattering lengths") {
        CHECK_THAT(protonated::get(form_factor_t::C).I0(),   WithinAbs(6.646, 1e-12));
        CHECK_THAT(protonated::get(form_factor_t::H).I0(),   WithinAbs(-3.739, 1e-12));
        CHECK_THAT(protonated::get(form_factor_t::CH2).I0(), WithinAbs(6.646 - 2*3.739, 1e-12));
        CHECK_THAT(deuterated::get(form_factor_t::H).I0(),   WithinAbs(6.671, 1e-12));
        CHECK_THAT(deuterated::get(form_factor_t::CH2).I0(), WithinAbs(6.646 + 2*6.671, 1e-12));
        CHECK_THAT(protonated::get(form_factor_t::NH3).I0(), WithinAbs(9.36 - 3*3.739, 1e-12));
    }

    SECTION("types without hydrogens are q-independent and isotope-independent") {
        for (auto type : {form_factor_t::C, form_factor_t::N, form_factor_t::O, form_factor_t::S, form_factor_t::OTHER}) {
            const auto& ffp = protonated::get(type);
            const auto& ffd = deuterated::get(type);
            CHECK_THAT(ffp.evaluate(1), WithinAbs(ffp.I0(), 1e-12));
            CHECK_THAT(ffd.evaluate(1), WithinAbs(ffp.I0(), 1e-12));
        }
    }

    SECTION("types with hydrogens vary with q") {
        for (auto type : {form_factor_t::CH, form_factor_t::CH2, form_factor_t::CH3, form_factor_t::NH, form_factor_t::NH2, form_factor_t::NH3, form_factor_t::OH, form_factor_t::SH}) {
            const auto& ff = protonated::get(type);
            CHECK(std::abs(ff.evaluate(1) - ff.I0()) > 0.1);
        }
    }

    SECTION("self-term equals the squared amplitude for single atoms") {
        for (auto type : {form_factor_t::H, form_factor_t::C, form_factor_t::N, form_factor_t::O, form_factor_t::S, form_factor_t::OTHER}) {
            const auto& ff = protonated::get(type);
            CHECK_THAT(ff.evaluate_self(1), WithinAbs(ff.evaluate(1)*ff.evaluate(1), 1e-12));
        }
    }

    SECTION("self-term is non-negative") {
        for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
            for (double q : {0., 0.25, 0.5, 1., 2., 5.}) {
                CHECK(protonated::get(static_cast<form_factor_t>(i)).evaluate_self(q) >= 0);
                CHECK(deuterated::get(static_cast<form_factor_t>(i)).evaluate_self(q) >= 0);
            }
        }
    }

    SECTION("excluded volume is not defined") {
        CHECK_THROWS(protonated::get(form_factor_t::EXCLUDED_VOLUME));
        CHECK_THROWS(deuterated::get(form_factor_t::EXCLUDED_VOLUME));
    }
}
