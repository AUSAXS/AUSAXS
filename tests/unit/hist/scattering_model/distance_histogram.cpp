#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distribution/Distribution1D.h>
#include <hist/scattering_model/DistanceHistogram.h>

#include <cmath>

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

TEST_CASE("DistanceHistogram::get_Rg") {
    SECTION("two point scatterers") {
        // two unit scatterers a distance d apart: p(0) = 2 (self terms), p(d) = 2. The transform carries an exp(-q^2) form
        // factor, i.e. a Gaussian of Rg^2 = 3 per scatterer, so Rg^2 = sum p r^2 / (2 sum p) + 3 = d^2/4 + 3
        int k = GENERATE(20, 80, 300);
        std::vector<double> counts(400, 0);
        counts[0] = 2;
        counts[k] = 2;
        hist::DistanceHistogram dhist{hist::Distribution1D(counts)};
        double d = dhist.get_d_axis()[k];
        CHECK_THAT(dhist.get_Rg(), Catch::Matchers::WithinAbs(std::sqrt(d*d/4 + 3), 1e-4));
    }

    SECTION("no contrast") {
        hist::DistanceHistogram dhist{hist::Distribution1D(std::vector<double>(100, 0))};
        CHECK(std::isnan(dhist.get_Rg()));
    }
}
