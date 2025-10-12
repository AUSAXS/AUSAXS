#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/decay/LinearDecay.h>
#include <rigidbody/parameters/decay/ExponentialDecay.h>

#include <cmath>

using namespace ausaxs;
using namespace rigidbody;

TEST_CASE("LinearDecay::factor") {
    auto val = GENERATE(10, 20, 30);
    double rate = 1./val;
    
    SECTION("Decay value " + std::to_string(val)) {
        parameter::decay::LinearDecay decay(val);
        for (unsigned int i = val; i > 0; --i) {
            REQUIRE_THAT(decay.next(), Catch::Matchers::WithinAbs(i*rate, 1e-6));
        }
    }
}

TEST_CASE("ExponentialDecay::factor") {
    auto val = GENERATE(10, 20, 30);
    
    SECTION("Decay value " + std::to_string(val)) {
        double rate = 2./val;
        parameter::decay::ExponentialDecay decay(val);
        for (double i = 0; i < val; ++i) {
            double factor = std::exp(-i*rate);
            REQUIRE_THAT(decay.next(), Catch::Matchers::WithinAbs(factor, 1e-6));
        }
    }
}