#include <catch2/catch_test_macros.hpp>

#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/distribution/Distribution1D.h>

using namespace ausaxs;

TEST_CASE("DistanceHistogram::get_weighted_counts") {
    SECTION("basic functionality") {
        std::vector<double> counts = {10, 20, 30, 40, 50};
        hist::Distribution1D dist(counts);
        hist::DistanceHistogram dhist(std::move(dist));
        
        const auto& weighted = dhist.get_weighted_counts();
        CHECK(weighted == counts);
        CHECK(weighted.size() == 5);
    }
}

TEST_CASE("DistanceHistogram::get_weighted_counts const-correctness") {
    SECTION("const histogram") {
        std::vector<double> counts = {10, 20, 30};
        hist::Distribution1D dist(counts);
        const hist::DistanceHistogram dhist(std::move(dist));
        
        const auto& weighted = dhist.get_weighted_counts();
        CHECK(weighted == counts);
    }

    SECTION("const reference") {
        std::vector<double> counts = {10, 20, 30};
        hist::Distribution1D dist(counts);
        hist::DistanceHistogram dhist_nonconst(std::move(dist));
        const hist::DistanceHistogram& dhist = dhist_nonconst;
        
        const auto& weighted = dhist.get_weighted_counts();
        CHECK(weighted == counts);
    }

    SECTION("returns const reference") {
        std::vector<double> counts = {10, 20, 30};
        hist::Distribution1D dist(counts);
        const hist::DistanceHistogram dhist(std::move(dist));
        
        const std::vector<double>& ref1 = dhist.get_weighted_counts();
        const std::vector<double>& ref2 = dhist.get_weighted_counts();
        CHECK(&ref1 == &ref2);
    }
}

TEST_CASE("DistanceHistogram::get_d_axis const-correctness") {
    SECTION("const histogram") {
        std::vector<double> counts = {10, 20, 30};
        hist::Distribution1D dist(counts);
        const hist::DistanceHistogram dhist(std::move(dist));
        
        const auto& axis = dhist.get_d_axis();
        CHECK(axis.size() == 3);
    }
}
