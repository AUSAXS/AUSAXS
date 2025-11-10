#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/xray/FormFactor.h>
#include <form_factor/xray/ExvFormFactor.h>
#include <form_factor/tables/FormFactorTableXray.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactor::constructor") {
    SECTION("5-Gaussian approximation") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        CHECK_THAT(ff.I0(), Catch::Matchers::WithinAbs(15.5, 1e-10));
    }

    SECTION("ExvFormFactor conversion") {
        ExvFormFactor exv(10.0);
        FormFactor ff(std::move(exv));
        CHECK(ff.I0() > 0);
    }
}

TEST_CASE("FormFactor::evaluate") {
    SECTION("at q = 0") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("at non-zero q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        double result = ff.evaluate(0.1);
        CHECK(result > 0);
        CHECK(result <= 1.0);
    }

    SECTION("decreases with q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        double val0 = ff.evaluate(0.0);
        double val1 = ff.evaluate(0.5);
        double val2 = ff.evaluate(1.0);
        
        CHECK(val0 >= val1);
        CHECK(val1 >= val2);
    }
}

TEST_CASE("FormFactor::I0") {
    SECTION("sum of coefficients") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        double expected = 1.0 + 2.0 + 3.0 + 4.0 + 5.0 + 0.5;
        CHECK_THAT(ff.I0(), Catch::Matchers::WithinAbs(expected, 1e-10));
    }

    SECTION("zero coefficients") {
        std::array<double, 5> a = {0, 0, 0, 0, 0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0;
        
        FormFactor ff(a, b, c);
        CHECK(ff.I0() == 0);
    }
}

TEST_CASE("FormFactor::set_normalization") {
    SECTION("change normalization") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
        
        double i0 = ff.I0();
        ff.set_normalization(2.0);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(i0 * 2.0, 1e-6));
    }

    SECTION("zero normalization") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        ff.set_normalization(0.0);
        CHECK(ff.evaluate(0) == 0.0);
    }
}

TEST_CASE("FormFactor::storage::atomic") {
    SECTION("get_form_factor for all types") {
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            const FormFactor& ff = storage::atomic::get_form_factor(static_cast<form_factor_t>(i));
            CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
            CHECK(ff.I0() > 0);
        }
    }

    SECTION("hydrogen") {
        const FormFactor& H = storage::atomic::get_form_factor(form_factor_t::H);
        CHECK_THAT(H.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("carbon") {
        const FormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        CHECK_THAT(C.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("nitrogen") {
        const FormFactor& N = storage::atomic::get_form_factor(form_factor_t::N);
        CHECK_THAT(N.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("oxygen") {
        const FormFactor& O = storage::atomic::get_form_factor(form_factor_t::O);
        CHECK_THAT(O.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("sulfur") {
        const FormFactor& S = storage::atomic::get_form_factor(form_factor_t::S);
        CHECK_THAT(S.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("excluded_volume") {
        const FormFactor& exv = storage::atomic::get_form_factor(form_factor_t::EXCLUDED_VOLUME);
        CHECK(exv.evaluate(0) > 0);
    }

    SECTION("invalid type throws") {
        CHECK_THROWS(storage::atomic::get_form_factor(form_factor_t::UNKNOWN));
    }
}

TEST_CASE("FormFactor::storage::atomic groups") {
    SECTION("CH_sp3") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::CH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("CH2_sp3") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::CH2);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("CH3_sp3") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::CH3);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("NH") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::NH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("NH2") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::NH2);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("NH3") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::NH3);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("OH") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::OH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("SH") {
        const FormFactor& ff = storage::atomic::get_form_factor(form_factor_t::SH);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }
}

TEST_CASE("FormFactor::consistency with q_axis") {
    SECTION("evaluate across q_axis") {
        const FormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        for (const double& q : constants::axes::q_vals) {
            double val = C.evaluate(q);
            CHECK(val > 0);
            CHECK(val <= 1.0);
        }
    }
}
