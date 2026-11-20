#include <catch2/catch_test_macros.hpp>

#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution2D.h>

using namespace ausaxs;

TEST_CASE("WeightedDistribution2D::WeightedDistribution2D") {
    SECTION("default constructor") {
        hist::WeightedDistribution2D dist;
        CHECK(dist.size_x() == 0);
        CHECK(dist.size_y() == 0);
    }

    SECTION("size constructor") {
        hist::WeightedDistribution2D dist(5, 10);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 10);
    }
}
