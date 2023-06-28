#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/ImageStack.h>

TEST_CASE("ImageStack::ImageStack") {
    SECTION("ExistingFile") {}
    SECTION("std::vector<Image>&") {}
    CHECK(false);
}

TEST_CASE("ImageStack::fit") {}
TEST_CASE("ImageStack::cutoff_scan") {}
TEST_CASE("ImageStack::cutoff_scan_fit") {}
TEST_CASE("ImageStack::get_fitted_water_factors") {}
TEST_CASE("ImageStack::get_fitted_water_factors_dataset") {}