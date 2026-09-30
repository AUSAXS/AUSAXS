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
        CHECK(dist.get_content(0, 0) == 1);
        CHECK(dist.get_content(1, 1) == 2);
        CHECK(dist.get_content(2, 2) == 3);
    }
}

TEST_CASE("Distribution2D::get_content") {
    hist::Distribution2D dist(5, 10);
    dist.index(0, 0) = 10;
    dist.index(1, 2) = 20;
    dist.index(3, 5) = 30;

    CHECK(dist.get_content(0, 0) == 10);
    CHECK(dist.get_content(1, 2) == 20);
    CHECK(dist.get_content(3, 5) == 30);
}

TEST_CASE("Distribution2D::add_index") {
    hist::Distribution2D dist(5, 10);
    dist.add_index(0, 0, 1);
    dist.add_index(1, 1, 2);
    dist.add_index(2, 2, 3);
    
    CHECK(dist.get_content(0, 0) == 1);
    CHECK(dist.get_content(1, 1) == 2);
    CHECK(dist.get_content(2, 2) == 3);
}

TEST_CASE("Distribution2D::increment_linear_index") {
    SECTION("basic increment") {
        hist::Distribution2D dist(3, 4);
        // Linear indexing for 2D: index = i * size_y + j
        // For (0,0): linear = 0
        // For (0,1): linear = 1
        // For (1,0): linear = 4
        dist.increment_linear_index(0);
        dist.increment_linear_index(1);
        dist.increment_linear_index(1);
        dist.increment_linear_index(4);
        dist.increment_linear_index(4);
        dist.increment_linear_index(4);
        
        CHECK(dist.linear_index(0) == 1);
        CHECK(dist.linear_index(1) == 2);
        CHECK(dist.linear_index(4) == 3);
    }

    SECTION("template parameter increment") {
        hist::Distribution2D dist(3, 4);
        dist.increment_linear_index<2>(0);
        dist.increment_linear_index<3>(1);
        dist.increment_linear_index<5>(2);
        
        CHECK(dist.linear_index(0) == 2);
        CHECK(dist.linear_index(1) == 3);
        CHECK(dist.linear_index(2) == 5);
    }
}
