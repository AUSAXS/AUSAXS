#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/Distribution3D.h>
#include <settings/HistogramSettings.h>

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

TEST_CASE("WeightedDistribution3D::get_weights") {
    SECTION("basic weighted axis") {
        hist::WeightedDistribution3D dist(10, 10, 10);
        dist.increment_index(0, 0, 0, 0.0f);
        dist.increment_index<2>(1, 1, 1, 1.25f);
        dist.increment_index<4>(2, 2, 2, 2.5f);

        auto weighted_bins = dist.get_weights();
        REQUIRE_THAT(weighted_bins[0], Catch::Matchers::WithinAbs(0, 1e-3));
        REQUIRE_THAT(weighted_bins[1], Catch::Matchers::WithinAbs(1.25, 1e-3));
        REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2.5, 1e-3));
    }

    SECTION("empty bins should return default axis values") {
        hist::WeightedDistribution3D dist(3, 3, 5);
        auto weighted_bins = dist.get_weights();
        
        // For empty bins (count=0), should return the default bin center
        for (size_t i = 0; i < weighted_bins.size(); ++i) {
            double expected = (i + 0.5) * settings::axes::bin_width;
            REQUIRE_THAT(weighted_bins[i], Catch::Matchers::WithinAbs(expected, 1e-6));
        }
    }

    SECTION("multiple xy entries should average correctly") {
        hist::WeightedDistribution3D dist(2, 2, 10);
        
        // Add values to bin z=2: distances at different xy positions
        dist.increment_index(0, 0, 2, 2.0f);
        dist.increment_index(0, 1, 2, 2.5f);
        dist.increment_index(1, 0, 2, 3.0f);
        dist.increment_index(1, 1, 2, 3.5f);
        
        auto weighted_bins = dist.get_weights();
        // Average: (2.0 + 2.5 + 3.0 + 3.5) / 4 = 2.75
        REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2.75, 1e-3));
    }
}
