#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/Distribution1D.h>

using namespace ausaxs;

TEST_CASE("WeightedDistribution1D::WeightedDistribution1D") {
    SECTION("default constructor") {
        hist::WeightedDistribution1D dist;
        CHECK(dist.size() == 0);
    }

    SECTION("size constructor") {
        hist::WeightedDistribution1D dist(10);
        CHECK(dist.size() == 10);
    }

    SECTION("from Distribution1D") {
        std::vector<double> data{1, 2, 3, 4, 5};
        hist::Distribution1D d1(data);
        hist::WeightedDistribution1D dist(d1);
        CHECK(dist.size() == 5);
        CHECK(dist.get_content(0) == 1);
        CHECK(dist.get_content(1) == 2);
        CHECK(dist.get_content(4) == 5);
    }

    SECTION("from vector") {
        std::vector<double> data{1, 2, 3, 4, 5};
        hist::WeightedDistribution1D dist(data);
        CHECK(dist.size() == 5);
        CHECK(dist.get_content(0) == 1);
        CHECK(dist.get_content(1) == 2);
        CHECK(dist.get_content(4) == 5);
    }
}

TEST_CASE("WeightedDistribution1D::get_content") {
    SECTION("vector") {
        hist::WeightedDistribution1D dist(5);
        dist.set_content(0, 1);
        dist.set_content(1, 2);
        dist.set_content(2, 3);
        
        auto content = dist.get_content();
        CHECK(content[0] == 1);
        CHECK(content[1] == 2);
        CHECK(content[2] == 3);
    }

    SECTION("indexed access") {
        hist::WeightedDistribution1D dist(5);
        dist.set_content(0, 1);
        dist.set_content(1, 2);
        dist.set_content(2, 3);
        
        CHECK(dist.get_content(0) == 1);
        CHECK(dist.get_content(1) == 2);
        CHECK(dist.get_content(2) == 3);
    }
}

TEST_CASE("WeightedDistribution1D::set_content") {
    hist::WeightedDistribution1D dist(5);
    dist.set_content(0, 10);
    dist.set_content(1, 20);
    dist.set_content(2, 30);
    CHECK(dist.get_content(0) == 10);
    CHECK(dist.get_content(1) == 20);
    CHECK(dist.get_content(2) == 30);
}

TEST_CASE("WeightedDistribution1D::add") {
    hist::WeightedDistribution1D dist(10);
    auto width = constants::axes::d_axis.width();
    
    dist.add(0, 1);
    dist.add(width/2, 2);
    
    CHECK(dist.get_content(0) == 1);
    CHECK(dist.get_content(1) == 2);
}

TEST_CASE("WeightedDistribution1D::add_index") {
    hist::WeightedDistribution1D dist(10);
    hist::detail::WeightedEntry entry1(5, 1, 1.5);
    hist::detail::WeightedEntry entry2(10, 2, 3.0);
    
    dist.add_index(0, entry1);
    dist.add_index(1, entry2);
    
    CHECK(dist.get_content(0) == 5);
    CHECK(dist.get_content(1) == 10);
}

TEST_CASE("WeightedDistribution1D::clear") {
    hist::WeightedDistribution1D dist(5);
    dist.set_content(2, 10);
    CHECK(dist.get_content(2) == 10);
    
    dist.clear(2);
    CHECK(dist.get_content(2) == 0);
}

TEST_CASE("WeightedDistribution1D::as_vector") {
    hist::WeightedDistribution1D dist(5);
    dist.set_content(0, 1);
    dist.set_content(1, 2);
    dist.set_content(2, 3);
    
    auto vec = dist.as_vector();
    CHECK(vec[0] == 1);
    CHECK(vec[1] == 2);
    CHECK(vec[2] == 3);
}

TEST_CASE("WeightedDistribution1D::get_weighted_axis") {
    hist::WeightedDistribution1D dist(10);
    auto width = constants::axes::d_axis.width();
    dist.add(0, 1);
    dist.add(width/2, 2);

    auto weighted_bins = dist.get_weighted_axis();
    REQUIRE_THAT(weighted_bins[0], Catch::Matchers::WithinAbs(0, 1e-3));
    REQUIRE_THAT(weighted_bins[1], Catch::Matchers::WithinAbs(width/2, 1e-3));
    REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2*width, 1e-3));
}

TEST_CASE("WeightedDistribution1D::set_bin_centers") {
    hist::WeightedDistribution1D dist(5);
    std::vector<double> centers{1.5, 2.5, 3.5, 4.5, 5.5};
    dist.set_bin_centers(centers);
    
    auto weighted_axis = dist.get_weighted_axis();
    CHECK(weighted_axis.size() == 5);
}

TEST_CASE("WeightedDistribution1D::operator+=") {
    hist::WeightedDistribution1D dist1(5);
    hist::WeightedDistribution1D dist2(5);
    
    dist1.set_content(0, 1);
    dist1.set_content(1, 2);
    dist2.set_content(0, 3);
    dist2.set_content(1, 4);
    
    dist1 += dist2;
    CHECK(dist1.get_content(0) == 4);
    CHECK(dist1.get_content(1) == 6);
}

TEST_CASE("WeightedDistribution1D::operator-=") {
    hist::WeightedDistribution1D dist1(5);
    hist::WeightedDistribution1D dist2(5);
    
    dist1.set_content(0, 5);
    dist1.set_content(1, 6);
    dist2.set_content(0, 3);
    dist2.set_content(1, 4);
    
    dist1 -= dist2;
    CHECK(dist1.get_content(0) == 2);
    CHECK(dist1.get_content(1) == 2);
}

TEST_CASE("WeightedDistribution1D::operator*") {
    hist::WeightedDistribution1D dist(5);
    dist.set_content(0, 1);
    dist.set_content(1, 2);
    dist.set_content(2, 3);
    
    auto result = 2.0 * dist;
    CHECK(result.get_content(0) == 2);
    CHECK(result.get_content(1) == 4);
    CHECK(result.get_content(2) == 6);
}
