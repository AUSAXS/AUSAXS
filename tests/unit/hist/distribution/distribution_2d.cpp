#include <catch2/catch_test_macros.hpp>

#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution2D.h>

using namespace ausaxs;

TEST_CASE("Distribution2D::Distribution2D") {
    SECTION("default constructor") {
        hist::Distribution2D dist;
        CHECK(dist.size_x() == 0);
        CHECK(dist.size_y() == 0);
    }

    SECTION("size constructor") {
        hist::Distribution2D dist(10, 10);
        CHECK(dist.size_x() == 10);
        CHECK(dist.size_y() == 10);
    }

    SECTION("from WeightedDistribution2D") {
        hist::WeightedDistribution2D wdist(5, 10);
        wdist.index(0, 0).value = 1;
        wdist.index(1, 1).value = 2;
        wdist.index(2, 2).value = 3;

        hist::Distribution2D dist(wdist);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 10);
        CHECK(dist.index(0, 0) == 1);
        CHECK(dist.index(1, 1) == 2);
        CHECK(dist.index(2, 2) == 3);
    }
}
