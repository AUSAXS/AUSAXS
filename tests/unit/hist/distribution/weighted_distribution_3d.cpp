#include <catch2/catch_test_macros.hpp>

#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution3D.h>

using namespace ausaxs;

TEST_CASE("WeightedDistribution3D::WeightedDistribution3D") {
    SECTION("default constructor") {
        hist::WeightedDistribution3D dist;
        CHECK(dist.size_x() == 0);
        CHECK(dist.size_y() == 0);
        CHECK(dist.size_z() == 0);
    }

    SECTION("size constructor") {
        hist::WeightedDistribution3D dist(5, 5, 10);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 5);
        CHECK(dist.size_z() == 10);
    }
}
