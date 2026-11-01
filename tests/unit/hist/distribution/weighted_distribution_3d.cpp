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

    SECTION("from Distribution3D") {
        hist::Distribution3D d3(5, 5, 10);
        d3.index(0, 0, 0) = 1;
        d3.index(1, 1, 1) = 2;
        d3.index(2, 2, 2) = 3;

        hist::WeightedDistribution3D dist(d3);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 5);
        CHECK(dist.size_z() == 10);
        CHECK(dist.index(0, 0, 0).value == 1);
        CHECK(dist.index(1, 1, 1).value == 2);
        CHECK(dist.index(2, 2, 2).value == 3);
    }
}

TEST_CASE("WeightedDistribution3D::increment_bin") {
    hist::WeightedDistribution3D dist(5, 5, 10);
    
    dist.increment_index(0, 0, 0, 0.5f);
    dist.increment_index(1, 1, 1, 1.5f);
    
    CHECK(dist.index(0, 0, 0).value == 1);
    CHECK(dist.index(1, 1, 1).value == 1);
}
