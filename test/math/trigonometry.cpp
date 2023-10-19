#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math/Trigonometry.h>

TEST_CASE("math::exp") {
    for (double x = -700; x < 10; x += 0.1) {
        CHECK_THAT(math::exp<9>(x), Catch::Matchers::WithinAbs(std::exp(x), 1e-3));
    }
}