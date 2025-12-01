#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/ExvTable.h>
#include <constants/Constants.h>

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
        auto set = detail::ExvFormFactorSet(constants::exv::standard);
        CHECK(set.C.is_initialized());
        CHECK(set.N.is_initialized());
        CHECK(set.O.is_initialized());
        CHECK(set.S.is_initialized());
    }

    SECTION("from Traube set") {
        auto set = detail::ExvFormFactorSet(constants::exv::Traube);
        CHECK(set.H.is_initialized());
        CHECK(set.C.is_initialized());
    }

    SECTION("from vdw set") {
        auto set = detail::ExvFormFactorSet(constants::exv::vdw);
        CHECK(set.H.is_initialized());
        CHECK(set.C.is_initialized());
    }
}

TEST_CASE("ExvFormFactorSet::get") {
    SECTION("all standard types") {
        auto set = detail::ExvFormFactorSet(constants::exv::standard);
        
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
        auto set = detail::ExvFormFactorSet(constants::exv::standard);
        CHECK_THROWS(set.get(form_factor_t::EXCLUDED_VOLUME));
    }
}

TEST_CASE("lookup::exv::standard") {
    SECTION("standard set is accessible") {
        const auto& set = lookup::exv::standard;
        CHECK(set.C.is_initialized());
        CHECK(set.N.is_initialized());
        CHECK(set.O.is_initialized());
        CHECK(set.S.is_initialized());
    }

    SECTION("all form factors evaluate properly") {
        const auto& set = lookup::exv::standard;
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            const ExvFormFactor& exv = set.get(static_cast<form_factor_t>(i));
            if (exv.is_initialized()) {
                CHECK_THAT(exv.evaluate_normalized(0), Catch::Matchers::WithinAbs(1.0, 1e-10));
            }
        }
    }
}

TEST_CASE("constants::exv::ExvSet") {
    SECTION("Traube set values") {
        CHECK(constants::exv::Traube.H > 0);
        CHECK(constants::exv::Traube.C > 0);
        CHECK(constants::exv::Traube.N > 0);
        CHECK(constants::exv::Traube.O > 0);
        CHECK(constants::exv::Traube.S > 0);
    }

    SECTION("vdw set values") {
        CHECK(constants::exv::vdw.H > 0);
        CHECK(constants::exv::vdw.C > 0);
        CHECK(constants::exv::vdw.N > 0);
        CHECK(constants::exv::vdw.O > 0);
        CHECK(constants::exv::vdw.S > 0);
    }

    SECTION("Voronoi sets") {
        CHECK(constants::exv::Voronoi_implicit_H.C > 0);
        CHECK(constants::exv::Voronoi_explicit_H.H > 0);
    }

    SECTION("MinimumFluctuation sets") {
        CHECK(constants::exv::MinimumFluctuation_implicit_H.C > 0);
        CHECK(constants::exv::MinimumFluctuation_explicit_H.H >= 0);
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
