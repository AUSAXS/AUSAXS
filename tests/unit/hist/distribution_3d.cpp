#include <catch2/catch_test_macros.hpp>

#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution3D.h>

using namespace ausaxs;

TEST_CASE("Distribution3D::Distribution3D") {
    SECTION("default constructor") {
        hist::Distribution3D dist;
        CHECK(dist.size_x() == 0);
        CHECK(dist.size_y() == 0);
        CHECK(dist.size_z() == 0);
    }

    SECTION("size constructor") {
        hist::Distribution3D dist(5, 10, 15);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 10);
        CHECK(dist.size_z() == 15);
    }

    SECTION("from WeightedDistribution3D") {
        hist::WeightedDistribution3D wdist(5, 5, 10);
        wdist.index(0, 0, 0).value = 1;
        wdist.index(1, 1, 1).value = 2;
        wdist.index(2, 2, 2).value = 3;

        hist::Distribution3D dist(wdist);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 5);
        CHECK(dist.size_z() == 10);
        CHECK(dist.index(0, 0, 0) == 1);
        CHECK(dist.index(1, 1, 1) == 2);
        CHECK(dist.index(2, 2, 2) == 3);
    }
}

TEST_CASE("Distribution3D::add_index") {
    hist::Distribution3D dist(5, 5, 10);
    dist.add_index(0, 0, 0, 1);
    dist.add_index(1, 1, 1, 2);
    dist.add_index(2, 2, 2, 3);
    
    CHECK(dist.index(0, 0, 0) == 1);
    CHECK(dist.index(1, 1, 1) == 2);
    CHECK(dist.index(2, 2, 2) == 3);
}
