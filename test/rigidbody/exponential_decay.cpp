#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/decay/ExponentialDecay.h>

#include <cmath>

using namespace rigidbody;

TEST_CASE("ExponentialDecay::factor") {
    auto val = GENERATE(10, 20, 30);
    
    SECTION("Decay value " + std::to_string(val)) {
        double rate = 2./val;
        parameters::decay::ExponentialDecay decay(val);
        for (double i = 0; i < val; ++i) {
            double factor = std::exp(-i*rate);
            REQUIRE_THAT(decay.get_factor(), Catch::Matchers::WithinAbs(factor, 1e-6));
        }
    }
}