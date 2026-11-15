#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/neutron/NeutronFormFactor.h>
#include <form_factor/neutron/NeutronSolvent.h>
#include <settings/ScatteringSettings.h>

#include <cmath>

using namespace ausaxs;
using namespace form_factor;
using Catch::Matchers::WithinAbs;

TEST_CASE("NeutronFormFactor::get") {
    using namespace neutron;
    settings::scattering::neutron::solvent_density = solvent_density(0);

    SECTION("types without hydrogens are q-independent") {
        for (auto type : {form_factor_t::C, form_factor_t::N, form_factor_t::O, form_factor_t::S, form_factor_t::OTHER}) {
            const auto& ff = get(type);
            CHECK_THAT(ff.evaluate(1), WithinAbs(ff.I0(), 1e-12));
        }
    }

    SECTION("types with hydrogens vary with q") {
        for (auto type : {form_factor_t::CH, form_factor_t::CH2, form_factor_t::CH3, form_factor_t::NH, form_factor_t::NH2, form_factor_t::NH3, form_factor_t::OH, form_factor_t::SH, form_factor_t::WATER}) {
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
        CHECK_THAT(get(form_factor_t::CH3).evaluate(1), WithinAbs(scattering_length::C + 3*scattering_length::H*sinc(group::bond::CH3), 1e-12));
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

TEST_CASE("neutron::d2o_fraction") {
    SECTION("inverts solvent_density") {
        for (double x : {0., 0.25, 0.42, 1.}) {
            CHECK_THAT(neutron::d2o_fraction(neutron::solvent_density(x)), WithinAbs(x, 1e-12));
        }
    }

    SECTION("extrapolates outside the H2O/D2O range") {
        double step = neutron::solvent_density(1) - neutron::solvent_density(0);
        CHECK_THAT(neutron::d2o_fraction(neutron::solvent_density(0) - 0.5*step), WithinAbs(-0.5, 1e-12));
        CHECK_THAT(neutron::d2o_fraction(neutron::solvent_density(1) + 0.5*step), WithinAbs(1.5, 1e-12));
    }
}

TEST_CASE("neutron::deuterium_fraction") {
    using namespace neutron;
    constexpr double x = 0.7, efficiency = 0.6;

    SECTION("hydration water follows the solvent") {
        CHECK(deuterium_fraction(form_factor_t::WATER, x, efficiency) == x);
    }

    SECTION("hydrogens on N, O and S exchange") {
        for (auto type : {form_factor_t::NH, form_factor_t::NH2, form_factor_t::NH3, form_factor_t::OH, form_factor_t::SH}) {
            CHECK_THAT(deuterium_fraction(type, x, efficiency), WithinAbs(efficiency*x, 1e-12));
        }
    }

    SECTION("hydrogens on C never exchange") {
        for (auto type : {form_factor_t::CH, form_factor_t::CH2, form_factor_t::CH3}) {
            CHECK(deuterium_fraction(type, x, efficiency) == 0);
        }
    }
}

TEST_CASE("NeutronFormFactor::get: hydrogen isotopes") {
    using namespace neutron;
    namespace b = scattering_length;
    auto mixed = [] (double x) {return (1-x)*b::H + x*b::D;};

    SECTION("the constexpr core mixes the implicit hydrogens") {
        constexpr double x = 0.3, efficiency = 0.8;
        static_assert(get(form_factor_t::CH3, 1, 1).I0() == b::C + 3*b::H);
        CHECK_THAT(get(form_factor_t::WATER, x, efficiency).I0(), WithinAbs(b::O + 2*mixed(x), 1e-12));
        CHECK_THAT(get(form_factor_t::NH2, x, efficiency).I0(), WithinAbs(b::N + 2*mixed(efficiency*x), 1e-12));
        CHECK_THAT(get(form_factor_t::CH2, x, efficiency).I0(), WithinAbs(b::C + 2*b::H, 1e-12));
    }

    SECTION("the isotope composition follows the solvent density") {
        settings::scattering::neutron::exchange_efficiency = 1;
        settings::scattering::neutron::solvent_density = solvent_density(1);
        CHECK_THAT(get(form_factor_t::WATER).I0(), WithinAbs(b::O + 2*b::D, 1e-12));
        CHECK_THAT(get(form_factor_t::OH).I0(), WithinAbs(b::O + b::D, 1e-12));
        CHECK_THAT(get(form_factor_t::CH).I0(), WithinAbs(b::C + b::H, 1e-12));

        settings::scattering::neutron::solvent_density = solvent_density(0.4);
        CHECK_THAT(get(form_factor_t::WATER).I0(), WithinAbs(b::O + 2*mixed(0.4), 1e-12));
        CHECK_THAT(get(form_factor_t::NH3).I0(), WithinAbs(b::N + 3*mixed(0.4), 1e-12));

        settings::scattering::neutron::exchange_efficiency = 0.5;
        CHECK_THAT(get(form_factor_t::WATER).I0(), WithinAbs(b::O + 2*mixed(0.4), 1e-12));
        CHECK_THAT(get(form_factor_t::NH3).I0(), WithinAbs(b::N + 3*mixed(0.2), 1e-12));
    }

    SECTION("explicit hydrogens are never exchanged") {
        settings::scattering::neutron::exchange_efficiency = 1;
        settings::scattering::neutron::solvent_density = solvent_density(1);
        CHECK(get(form_factor_t::H).I0() == b::H);
        CHECK(get(form_factor_t::D).I0() == b::D);
    }
}
