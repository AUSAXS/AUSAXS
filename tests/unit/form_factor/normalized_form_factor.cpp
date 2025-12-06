#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/ExvFormFactor.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("NormalizedFormFactor::constructor") {
    SECTION("5-Gaussian approximation") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        NormalizedFormFactor ff(a, b, c);
        CHECK_THAT(ff.I0(), Catch::Matchers::WithinAbs(15.5, 1e-10));
    }

    SECTION("ExvFormFactor conversion") {
        ExvFormFactor exv(10.0);
        NormalizedFormFactor ff(std::move(exv));
        CHECK(ff.I0() > 0);
    }
}

TEST_CASE("NormalizedFormFactor::evaluate") {
    SECTION("at q = 0") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        NormalizedFormFactor ff(a, b, c);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("at non-zero q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        NormalizedFormFactor ff(a, b, c);
        double result = ff.evaluate(0.1);
        CHECK(result > 0);
        CHECK(result <= 1.0);
    }

    SECTION("decreases with q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        NormalizedFormFactor ff(a, b, c);
        double val0 = ff.evaluate(0.0);
        double val1 = ff.evaluate(0.5);
        double val2 = ff.evaluate(1.0);
        
        CHECK(val0 >= val1);
        CHECK(val1 >= val2);
    }
}

TEST_CASE("NormalizedFormFactor::I0") {
    SECTION("sum of coefficients") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        NormalizedFormFactor ff(a, b, c);
        double expected = 1.0 + 2.0 + 3.0 + 4.0 + 5.0 + 0.5;
        CHECK_THAT(ff.I0(), Catch::Matchers::WithinAbs(expected, 1e-10));
    }

    SECTION("zero coefficients") {
        std::array<double, 5> a = {0, 0, 0, 0, 0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0;
        
        NormalizedFormFactor ff(a, b, c);
        CHECK(ff.I0() == 0);
    }
}

TEST_CASE("NormalizedFormFactor::set_normalization") {
    SECTION("change normalization") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        NormalizedFormFactor ff(a, b, c);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
        
        double i0 = ff.I0();
        ff.set_normalization(2.0);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(i0 * 2.0, 1e-6));
    }

    SECTION("zero normalization") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        NormalizedFormFactor ff(a, b, c);
        ff.set_normalization(0.0);
        CHECK(ff.evaluate(0) == 0.0);
    }
}

TEST_CASE("NormalizedFormFactor::storage::atomic") {
    SECTION("get_form_factor for all types") {
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            const NormalizedFormFactor& ff = lookup::atomic::normalized::get(static_cast<form_factor_t>(i));
            CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
            CHECK(ff.I0() > 0);
        }
    }

    SECTION("hydrogen") {
        const NormalizedFormFactor& H = lookup::atomic::normalized::get(form_factor_t::H);
        CHECK_THAT(H.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("carbon") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        CHECK_THAT(C.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("nitrogen") {
        const NormalizedFormFactor& N = lookup::atomic::normalized::get(form_factor_t::N);
        CHECK_THAT(N.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("oxygen") {
        const NormalizedFormFactor& O = lookup::atomic::normalized::get(form_factor_t::O);
        CHECK_THAT(O.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("sulfur") {
        const NormalizedFormFactor& S = lookup::atomic::normalized::get(form_factor_t::S);
        CHECK_THAT(S.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("excluded_volume") {
        const NormalizedFormFactor& exv = lookup::atomic::normalized::get(form_factor_t::EXCLUDED_VOLUME);
        CHECK(exv.evaluate(0) > 0);
    }

    SECTION("invalid type throws") {
        CHECK_THROWS(lookup::atomic::normalized::get(form_factor_t::UNKNOWN));
    }
}

TEST_CASE("NormalizedFormFactor::storage::atomic groups") {
    SECTION("CH_sp3") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::CH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("CH2_sp3") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::CH2);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("CH3_sp3") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::CH3);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("NH") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::NH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("NH2") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::NH2);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("NH3") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::NH3);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("OH") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::OH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("SH") {
        const NormalizedFormFactor& ff = lookup::atomic::normalized::get(form_factor_t::SH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }
}

TEST_CASE("NormalizedFormFactor::consistency with q_axis") {
    SECTION("evaluate across q_axis") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        for (const double& q : constants::axes::q_vals) {
            double val = C.evaluate(q);
            CHECK(val > 0);
            CHECK(val <= 1.0);
        }
    }
}
