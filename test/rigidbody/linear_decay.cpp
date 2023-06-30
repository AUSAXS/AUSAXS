#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/decay/LinearDecay.h>

using namespace rigidbody;

TEST_CASE("LinearDecay::factor") {
    auto val = GENERATE(10, 20, 30);
    double rate = 1./val;
    
    SECTION("Decay value " + std::to_string(val)) {
        parameters::decay::LinearDecay decay(val);
        for (unsigned int i = val; i > 0; --i) {
            REQUIRE_THAT(decay.get_factor(), Catch::Matchers::WithinAbs(i*rate, 1e-6));
        }
    }
}