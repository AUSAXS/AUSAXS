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

TEST_CASE("Distribution2D::add") {
    hist::Distribution2D dist(5, 10);
    dist.add(0, 0, 1);
    dist.add(1, 1.5, 2);
    dist.add(2, 3.3, 3);
    
    CHECK(dist.get_content(0, 0) == 1);
    CHECK(dist.get_content(1, 2) == 2);
    CHECK(dist.get_content(2, 3) == 3);
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
