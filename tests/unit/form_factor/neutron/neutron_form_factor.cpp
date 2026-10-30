#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/neutron/NeutronFormFactor.h>

#include <cmath>

using namespace ausaxs;
using namespace form_factor;
using Catch::Matchers::WithinAbs;

TEST_CASE("NeutronFormFactor::get") {
    using namespace neutron;

    SECTION("q = 0 values are the summed scattering lengths") {
        CHECK_THAT(get(form_factor_t::C).I0(),   WithinAbs(6.646, 1e-12));
        CHECK_THAT(get(form_factor_t::H).I0(),   WithinAbs(-3.739, 1e-12));
        CHECK_THAT(get(form_factor_t::CH2).I0(), WithinAbs(6.646 - 2*3.739, 1e-12));
        CHECK_THAT(get(form_factor_t::NH3).I0(), WithinAbs(9.36 - 3*3.739, 1e-12));
    }

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

    SECTION("groups place their hydrogens at the nuclear X-H distance") {
        auto sinc = [] (double x) {return std::sin(x)/x;};
        CHECK_THAT(get(form_factor_t::CH3).evaluate(1), WithinAbs(6.646 - 3*3.739*sinc(group::bond::CH), 1e-12));
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
