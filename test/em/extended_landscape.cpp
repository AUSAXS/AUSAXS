#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/detail/ExtendedLandscape.h>

TEST_CASE("ExtendedLandscape::ExtendedLandscape") {
    mini::Evaluation eval1({1, 2, 3}, 4);
    mini::Evaluation eval2({5, 6, 7}, 8);
    mini::Landscape landscape1({eval1, eval2});

    SECTION("double, double, Landscape&&") {
        em::detail::ExtendedLandscape landscape(1, 2, std::move(landscape1));
        CHECK(landscape.cutoff == 1);
        CHECK(landscape.mass == 2);
        REQUIRE(landscape.strip.evals.size() == 2);
        CHECK(landscape.strip.evals[0] == eval1);
        CHECK(landscape.strip.evals[1] == eval2);
    }
}