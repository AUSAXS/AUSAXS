#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/xray/XrayFormFactor.h>

#include <cmath>

using namespace ausaxs;
using namespace form_factor;
using Catch::Matchers::WithinAbs;

TEST_CASE("AtomicFormFactor::constructor") {
    SECTION("5-Gaussian approximation") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        xray::AtomicFormFactor ff(a, b, c);
        CHECK_THAT(ff.I0(), Catch::Matchers::WithinAbs(15.5, 1e-10));
    }

    SECTION("ExvFormFactor conversion") {
        auto exv = factory::construct_exv_form_factor(10.0);
        xray::AtomicFormFactor ff(std::move(exv));
        CHECK(ff.I0() > 0);
    }
}

TEST_CASE("AtomicFormFactor::evaluate") {
    SECTION("at q = 0") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        xray::AtomicFormFactor ff(a, b, c);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(15.5, 1e-6));
    }

    SECTION("at non-zero q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        xray::AtomicFormFactor ff(a, b, c);
        double result = ff.evaluate(0.1);
        CHECK(result > 0);
        CHECK(result < 15.5);
    }

    SECTION("decreases with q") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        xray::AtomicFormFactor ff(a, b, c);
        double q1 = ff.evaluate(0.1);
        double q2 = ff.evaluate(0.5);
        double q3 = ff.evaluate(1.0);
        
        CHECK(q1 > q2);
        CHECK(q2 > q3);
    }
}

TEST_CASE("AtomicFormFactor::set_normalization") {
    SECTION("manual normalization") {
        std::array<double, 5> a = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::array<double, 5> b = {0.1, 0.2, 0.3, 0.4, 0.5};
        double c = 0.5;
        
        xray::AtomicFormFactor ff(a, b, c);
        ff.set_normalization(2);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(2, 1e-9));

        ff.set_normalization(30);
        CHECK_THAT(ff.evaluate(0), Catch::Matchers::WithinAbs(30, 1e-9));
    }
}

TEST_CASE("FormFactor::xray::raw") {
    SECTION("get hydrogen form factor") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::H);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get carbon form factor") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::C);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get nitrogen form factor") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::N);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get oxygen form factor") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::O);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get sulfur form factor") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::S);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get excluded volume form factor") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::EXCLUDED_VOLUME);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("get form factor by type") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::C);
        CHECK(ff.I0() > 0);
        CHECK(ff.evaluate(0) > 0);
    }
}

TEST_CASE("FormFactor::comparison with normalized") {
    SECTION("raw form factors are not normalized to 1") {
        const xray::FormFactor& ff_h = xray::raw::get(form_factor_t::H);
        const xray::FormFactor& ff_c = xray::raw::get(form_factor_t::C);
        const xray::FormFactor& ff_n = xray::raw::get(form_factor_t::N);
        const xray::FormFactor& ff_o = xray::raw::get(form_factor_t::O);
        
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

TEST_CASE("FormFactor::xray::groups") {
    using namespace xray;
    auto sinc = [] (double x) {return std::sin(x)/x;};

    SECTION("q = 0 values are the summed electron counts") {
        for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
            auto type = static_cast<form_factor_t>(i);
            INFO("form factor " << to_string(type));
            CHECK_THAT(raw::get(type).I0(), WithinAbs(constants::charge::nuclear::get_charge(type), 0.01));
        }
    }

    SECTION("groups are built from the atomic form factors of their atoms") {
        const AtomicFormFactor c(coefficients::C), h(coefficients::H);
        for (double q : {0.1, 0.5, 1., 2.}) {
            CHECK_THAT(raw::get(form_factor_t::CH3).evaluate(q), WithinAbs(c.evaluate(q) + 3*h.evaluate(q)*sinc(q*group::bond::CH), 1e-12));
        }
    }

    SECTION("self-term equals the squared amplitude for single atoms") {
        for (auto type : {form_factor_t::EXCLUDED_VOLUME, form_factor_t::H, form_factor_t::C, form_factor_t::N, form_factor_t::O, form_factor_t::S, form_factor_t::OTHER}) {
            const auto& ff = raw::get(type);
            CHECK_THAT(ff.evaluate_self(1), WithinAbs(ff.evaluate(1)*ff.evaluate(1), 1e-12));
        }
    }

    SECTION("self-term of a group exceeds the squared amplitude") {
        // the difference is the variance of the amplitude over the orientations of the hydrogens, which is zero only at q = 0
        for (auto type : {form_factor_t::CH, form_factor_t::CH2, form_factor_t::CH3, form_factor_t::NH, form_factor_t::NH2, form_factor_t::NH3, form_factor_t::OH, form_factor_t::SH}) {
            const auto& ff = raw::get(type);
            INFO("form factor " << to_string(type));
            CHECK_THAT(ff.evaluate_self(0), WithinAbs(ff.I0()*ff.I0(), 1e-9));
            for (double q : {0.5, 1., 2., 5.}) {
                CHECK(ff.evaluate(q)*ff.evaluate(q) < ff.evaluate_self(q));
            }
        }
    }
}
