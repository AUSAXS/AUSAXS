#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/ExvFormFactor.h>
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
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(15.5, 1e-6));
    }

    SECTION("at non-zero q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        double result = ff.evaluate(0.1);
        CHECK(result > 0);
        CHECK(result < 15.5);
    }

    SECTION("decreases with q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        double q1 = ff.evaluate(0.1);
        double q2 = ff.evaluate(0.5);
        double q3 = ff.evaluate(1.0);
        
        CHECK(q1 > q2);
        CHECK(q2 > q3);
    }
}

TEST_CASE("FormFactor::set_normalization") {
    SECTION("manual normalization") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        FormFactor ff(a, b, c);
        ff.set_normalization(2.0);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(31.0, 1e-6));
    }
}

TEST_CASE("FormFactor::lookup::atomic::raw") {
    SECTION("get hydrogen form factor") {
        const FormFactor& ff = lookup::atomic::raw::H;
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get carbon form factor") {
        const FormFactor& ff = lookup::atomic::raw::C;
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get nitrogen form factor") {
        const FormFactor& ff = lookup::atomic::raw::N;
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get oxygen form factor") {
        const FormFactor& ff = lookup::atomic::raw::O;
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get sulfur form factor") {
        const FormFactor& ff = lookup::atomic::raw::S;
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get excluded volume form factor") {
        const FormFactor& ff = lookup::atomic::raw::excluded_volume;
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get form factor by type") {
        const FormFactor& ff = lookup::atomic::raw::get(form_factor_t::C);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }
}

TEST_CASE("FormFactor::comparison with normalized") {
    SECTION("raw form factors are not normalized to 1") {
        const FormFactor& ff_h = lookup::atomic::raw::H;
        const FormFactor& ff_c = lookup::atomic::raw::C;
        const FormFactor& ff_n = lookup::atomic::raw::N;
        const FormFactor& ff_o = lookup::atomic::raw::O;
        
        CHECK(ff_h.evaluate(0) != 1.0);
        CHECK(ff_c.evaluate(0) != 1.0);
        CHECK(ff_n.evaluate(0) != 1.0);
        CHECK(ff_o.evaluate(0) != 1.0);
        
        CHECK(ff_h.evaluate(0) > 0);
        CHECK(ff_c.evaluate(0) > 0);
        CHECK(ff_n.evaluate(0) > 0);
        CHECK(ff_o.evaluate(0) > 0);
    }
}
