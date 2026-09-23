#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/ExvTable.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <numbers>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("ExvFormFactor::constructor") {
    SECTION("positive volume") {
        double volume = 10.0;
        ExvFormFactor exv(volume);
        CHECK(exv.is_initialized());
        CHECK(exv.exponent > 0);
        CHECK(exv.q0 > 0);
    }

    SECTION("zero volume") {
        double volume = 0.0;
        ExvFormFactor exv(volume);
        CHECK_FALSE(exv.is_initialized());
        CHECK(exv.exponent == 0);
    }

    SECTION("large volume") {
        double volume = 1000.0;
        ExvFormFactor exv(volume);
        CHECK(exv.is_initialized());
        CHECK(exv.exponent > 0);
        CHECK(exv.q0 > 0);
    }
}

TEST_CASE("ExvFormFactor::evaluate") {
    SECTION("at q = 0") {
        double volume = 10.0;
        ExvFormFactor exv(volume);
        double val = exv.evaluate(0);
        CHECK_THAT(val, Catch::Matchers::WithinAbs(exv.q0, 1e-10));
    }

    SECTION("decreases with q") {
        double volume = 10.0;
        ExvFormFactor exv(volume);
        double val0 = exv.evaluate(0.0);
        double val1 = exv.evaluate(0.5);
        double val2 = exv.evaluate(1.0);
        
        CHECK(val0 >= val1);
        CHECK(val1 >= val2);
    }

    SECTION("positive values") {
        double volume = 10.0;
        ExvFormFactor exv(volume);
        for (double q = 0; q < 2.0; q += 0.1) {
            CHECK(exv.evaluate(q) > 0);
        }
    }
}

TEST_CASE("ExvFormFactor::evaluate_normalized") {
    SECTION("at q = 0") {
        double volume = 10.0;
        ExvFormFactor exv(volume);
        CHECK_THAT(exv.evaluate_normalized(0), Catch::Matchers::WithinAbs(1.0, 1e-10));
    }

    SECTION("decreases with q") {
        double volume = 10.0;
        ExvFormFactor exv(volume);
        double val0 = exv.evaluate_normalized(0.0);
        double val1 = exv.evaluate_normalized(0.5);
        double val2 = exv.evaluate_normalized(1.0);
        
        CHECK_THAT(val0, Catch::Matchers::WithinAbs(1.0, 1e-10));
        CHECK(val1 < val0);
        CHECK(val2 < val1);
    }

    SECTION("less than or equal to 1") {
        double volume = 10.0;
        ExvFormFactor exv(volume);
        for (double q = 0; q < 2.0; q += 0.1) {
            CHECK(exv.evaluate_normalized(q) <= 1.0);
            CHECK(exv.evaluate_normalized(q) > 0);
        }
    }
}

TEST_CASE("ExvFormFactor::is_initialized") {
    SECTION("initialized") {
        ExvFormFactor exv(10.0);
        CHECK(exv.is_initialized());
    }

    SECTION("not initialized") {
        ExvFormFactor exv(0.0);
        CHECK_FALSE(exv.is_initialized());
    }
}

TEST_CASE("ExvFormFactor::default_values") {
    SECTION("zero volume gives uninitialized") {
        ExvFormFactor exv(0.0);
        CHECK(exv.exponent == 0);
        CHECK_FALSE(exv.is_initialized());
    }
}

TEST_CASE("ExvFormFactor::volume_relationship") {
    SECTION("larger volume means larger q0") {
        ExvFormFactor exv1(10.0);
        ExvFormFactor exv2(20.0);
        CHECK(exv2.q0 > exv1.q0);
    }

    SECTION("larger volume means larger exponent") {
        ExvFormFactor exv1(10.0);
        ExvFormFactor exv2(20.0);
        CHECK(exv2.exponent > exv1.exponent);
    }
}

TEST_CASE("ExvFormFactorSet::constructor") {
    SECTION("from standard set") {
        auto set = detail::ExvFormFactorSet(ExvTableManager::get_default_exv_table());
        CHECK(set.get(form_factor_t::C).is_initialized());
        CHECK(set.get(form_factor_t::N).is_initialized());
        CHECK(set.get(form_factor_t::O).is_initialized());
        CHECK(set.get(form_factor_t::S).is_initialized());
    }

    SECTION("from Traube set") {
        auto set = detail::ExvFormFactorSet(constants::exv::Traube);
        CHECK(set.get(form_factor_t::H).is_initialized());
        CHECK(set.get(form_factor_t::C).is_initialized());
    }

    SECTION("from vdw set") {
        auto set = detail::ExvFormFactorSet(constants::exv::vdw);
        CHECK(set.get(form_factor_t::H).is_initialized());
        CHECK(set.get(form_factor_t::C).is_initialized());
    }
}

TEST_CASE("ExvFormFactorSet::get") {
    SECTION("all standard types") {
        auto set = detail::ExvFormFactorSet(ExvTableManager::get_default_exv_table());
        
        CHECK(set.get(form_factor_t::C).is_initialized());
        CHECK(set.get(form_factor_t::CH).is_initialized());
        CHECK(set.get(form_factor_t::CH2).is_initialized());
        CHECK(set.get(form_factor_t::CH3).is_initialized());
        CHECK(set.get(form_factor_t::N).is_initialized());
        CHECK(set.get(form_factor_t::NH).is_initialized());
        CHECK(set.get(form_factor_t::NH2).is_initialized());
        CHECK(set.get(form_factor_t::NH3).is_initialized());
        CHECK(set.get(form_factor_t::O).is_initialized());
        CHECK(set.get(form_factor_t::OH).is_initialized());
        CHECK(set.get(form_factor_t::S).is_initialized());
        CHECK(set.get(form_factor_t::SH).is_initialized());
        CHECK(set.get(form_factor_t::OTHER).is_initialized());
    }

    SECTION("invalid type throws") {
        auto set = detail::ExvFormFactorSet(ExvTableManager::get_default_exv_table());
        CHECK_THROWS(set.get(form_factor_t::EXCLUDED_VOLUME));
    }
}

TEST_CASE("ExvTableManager::get_current_exv_form_factor_set") {
    SECTION("standard set is accessible") {
        const auto& set = ExvTableManager::get_current_exv_form_factor_set();
        CHECK(set.get(form_factor_t::C).is_initialized());
        CHECK(set.get(form_factor_t::N).is_initialized());
        CHECK(set.get(form_factor_t::O).is_initialized());
        CHECK(set.get(form_factor_t::S).is_initialized());
    }

    SECTION("all form factors evaluate properly") {
        const auto& set = ExvTableManager::get_current_exv_form_factor_set();
        for (int i = 1; i < total_ff_count; ++i) {
            const ExvFormFactor& exv = set.get(static_cast<form_factor_t>(i));
            if (exv.is_initialized()) {
                CHECK_THAT(exv.evaluate_normalized(0), Catch::Matchers::WithinAbs(1.0, 1e-10));
            }
        }
    }
}

TEST_CASE("constants::exv::ExvSet") {
    SECTION("Traube set values") {
        CHECK(constants::exv::Traube.get(form_factor_t::H) > 0);
        CHECK(constants::exv::Traube.get(form_factor_t::C) > 0);
        CHECK(constants::exv::Traube.get(form_factor_t::N) > 0);
        CHECK(constants::exv::Traube.get(form_factor_t::O) > 0);
        CHECK(constants::exv::Traube.get(form_factor_t::S) > 0);
    }

    SECTION("vdw set values") {
        CHECK(constants::exv::vdw.get(form_factor_t::H) > 0);
        CHECK(constants::exv::vdw.get(form_factor_t::C) > 0);
        CHECK(constants::exv::vdw.get(form_factor_t::N) > 0);
        CHECK(constants::exv::vdw.get(form_factor_t::O) > 0);
        CHECK(constants::exv::vdw.get(form_factor_t::S) > 0);
    }

    SECTION("Voronoi sets") {
        CHECK(constants::exv::Voronoi_implicit_H.get(form_factor_t::C) > 0);
        CHECK(constants::exv::Voronoi_explicit_H.get(form_factor_t::H) > 0);
    }

    SECTION("MinimumFluctuation sets") {
        CHECK(constants::exv::MinimumFluctuation_implicit_H.get(form_factor_t::C) > 0);
        CHECK(constants::exv::MinimumFluctuation_explicit_H.get(form_factor_t::H) >= 0);
    }
}

TEST_CASE("constants::exv::volume") {
    SECTION("sphere volume") {
        double radius = 1.0;
        double volume = constants::exv::detail::volume(radius);
        double expected = 4.0 * std::numbers::pi / 3.0;
        CHECK_THAT(volume, Catch::Matchers::WithinRel(expected, 1e-10));
    }

    SECTION("zero radius") {
        double volume = constants::exv::detail::volume(0.0);
        CHECK(volume == 0.0);
    }

    SECTION("larger radius") {
        double volume1 = constants::exv::detail::volume(1.0);
        double volume2 = constants::exv::detail::volume(2.0);
        CHECK(volume2 > volume1);
    }
}

TEST_CASE("constants::exv::detail::ExvSet: optional entries") {
    SECTION("built-in sets match the tabulated values") {
        CHECK_THAT(constants::exv::Traube.get(form_factor_t::CH3), Catch::Matchers::WithinRel(31.89, 1e-12));
        CHECK(constants::exv::Voronoi_implicit_H.get(form_factor_t::NH2) == 22.129);
        CHECK(constants::exv::MinimumFluctuation_explicit_H.get(form_factor_t::SH) == 28.475);
        CHECK(constants::exv::vdw.get(form_factor_t::OTHER) == constants::exv::Ar);
    }

    SECTION("the excluded volume type has no entry") {
        CHECK_FALSE(constants::exv::Traube.contains(form_factor_t::EXCLUDED_VOLUME));
        CHECK_THROWS(constants::exv::Traube.get(form_factor_t::EXCLUDED_VOLUME));
    }

    SECTION("missing entries propagate to the form factor set") {
        auto set = constants::exv::Traube;
        set.volumes[static_cast<int>(form_factor_t::NH)].reset();
        CHECK_FALSE(set.contains(form_factor_t::NH));
        CHECK_THROWS(set.get(form_factor_t::NH));

        auto ffset = detail::ExvFormFactorSet(set);
        CHECK_FALSE(ffset.contains(form_factor_t::NH));
        CHECK_THROWS(ffset.get(form_factor_t::NH));
        CHECK(ffset.contains(form_factor_t::N));
    }
}
