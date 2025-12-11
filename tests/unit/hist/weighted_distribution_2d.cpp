#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/Distribution2D.h>

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

    SECTION("from Distribution2D") {
        hist::Distribution2D d2(5, 10);
        d2.index(0, 0) = 1;
        d2.index(1, 1) = 2;
        d2.index(2, 2) = 3;

        hist::WeightedDistribution2D dist(d2);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 10);
        CHECK(dist.index(0, 0).value == 1);
        CHECK(dist.index(1, 1).value == 2);
        CHECK(dist.index(2, 2).value == 3);
    }
}

TEST_CASE("WeightedDistribution2D::increment_bin") {
    hist::WeightedDistribution2D dist(5, 10);
    
    dist.increment_bin(0, 0, 0.5f);
    dist.increment_bin(1, 1, 1.5f);
    
    CHECK(dist.index(0, 0).value == 1);
    CHECK(dist.index(1, 1).value == 1);
}

TEST_CASE("WeightedDistribution2D::get_weighted_axis") {
    hist::WeightedDistribution2D dist(10, 10);
    dist.increment_bin(0, 0, 0.0f);
    dist.increment_bin<2>(1, 1, 1.25f);
    dist.increment_bin<4>(2, 2, 2.5f);

    auto weighted_bins = dist.get_weighted_axis();
    REQUIRE_THAT(weighted_bins[0], Catch::Matchers::WithinAbs(0, 1e-3));
    REQUIRE_THAT(weighted_bins[1], Catch::Matchers::WithinAbs(1.25, 1e-3));
    REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2.5, 1e-3));
}
