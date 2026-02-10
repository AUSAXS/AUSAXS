#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/parameters/decay/LinearDecay.h>
#include <rigidbody/parameters/decay/ExponentialDecay.h>
#include <rigidbody/parameters/decay/NoDecay.h>

using namespace ausaxs::rigidbody::parameter::decay;

TEST_CASE("NoDecay::next") {
    NoDecay decay;

    SECTION("always returns 1") {
        for (int i = 0; i < 100; ++i) {
            REQUIRE(decay.next() == 1.0);
        }
    }
}

TEST_CASE("LinearDecay::next") {
    SECTION("starts at 1") {
        LinearDecay decay(100);
        REQUIRE_THAT(decay.next(), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("monotonically decreasing") {
        LinearDecay decay(100);
        double prev = decay.next();
        for (int i = 1; i < 50; ++i) {
            double cur = decay.next();
            REQUIRE(cur < prev);
            prev = cur;
        }
    }

    SECTION("set_iterations") {
        LinearDecay decay(100);
        decay.set_iterations(200);
        REQUIRE(decay.get_iterations() == 200);
        double first = decay.next();
        double second = decay.next();
        REQUIRE(second < first);
    }
}

TEST_CASE("ExponentialDecay::next") {
    SECTION("starts at 1") {
        ExponentialDecay decay(100);
        REQUIRE_THAT(decay.next(), Catch::Matchers::WithinAbs(1.0, 1e-6));
    }

    SECTION("monotonically decreasing") {
        ExponentialDecay decay(100);
        double prev = decay.next();
        for (int i = 1; i < 50; ++i) {
            double cur = decay.next();
            REQUIRE(cur < prev);
            prev = cur;
        }
    }

    SECTION("always positive") {
        ExponentialDecay decay(100);
        for (int i = 0; i < 200; ++i) {
            REQUIRE(decay.next() > 0);
        }
    }

    SECTION("set_iterations") {
        ExponentialDecay decay(100);
        decay.set_iterations(200);
        REQUIRE(decay.get_iterations() == 200);
    }
}
