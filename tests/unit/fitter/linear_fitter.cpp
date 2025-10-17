#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fitter/LinearFitter.h>
#include <dataset/SimpleDataset.h>
#include <hist/intensity_calculator/DistanceHistogram.h>

using namespace ausaxs;
using namespace ausaxs::fitter;

TEST_CASE("LinearFitter::constructor") {
    std::vector<double> x = {0.1, 0.2, 0.3, 0.4};
    std::vector<double> y = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> yerr = {0.1, 0.1, 0.1, 0.1};
    SimpleDataset data(x, y, yerr);
    
    LinearFitter fitter(data);
    REQUIRE(fitter.size() == 4);
    REQUIRE(fitter.dof() == 2);
}

TEST_CASE("LinearFitter::get_data") {
    std::vector<double> x = {0.1, 0.2, 0.3};
    std::vector<double> y = {1.0, 2.0, 3.0};
    std::vector<double> yerr = {0.1, 0.1, 0.1};
    SimpleDataset data(x, y, yerr);
    
    LinearFitter fitter(data);
    SimpleDataset retrieved = fitter.get_data();
    
    REQUIRE(retrieved.size() == 3);
    REQUIRE(retrieved == data);
}

TEST_CASE("LinearFitter::size") {
    SECTION("3 points") {
        std::vector<double> x = {0.1, 0.2, 0.3};
        std::vector<double> y = {1.0, 2.0, 3.0};
        std::vector<double> yerr = {0.1, 0.1, 0.1};
        SimpleDataset data(x, y, yerr);
        
        LinearFitter fitter(data);
        REQUIRE(fitter.size() == 3);
    }
    
    SECTION("5 points") {
        std::vector<double> x = {0.1, 0.2, 0.3, 0.4, 0.5};
        std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::vector<double> yerr = {0.1, 0.1, 0.1, 0.1, 0.1};
        SimpleDataset data(x, y, yerr);
        
        LinearFitter fitter(data);
        REQUIRE(fitter.size() == 5);
    }
}

TEST_CASE("LinearFitter::dof") {
    SECTION("3 points") {
        std::vector<double> x = {0.1, 0.2, 0.3};
        std::vector<double> y = {1.0, 2.0, 3.0};
        std::vector<double> yerr = {0.1, 0.1, 0.1};
        SimpleDataset data(x, y, yerr);
        
        LinearFitter fitter(data);
        REQUIRE(fitter.dof() == 1);
    }
    
    SECTION("10 points") {
        std::vector<double> x(10);
        std::vector<double> y(10);
        std::vector<double> yerr(10, 0.1);
        for (int i = 0; i < 10; ++i) {
            x[i] = 0.1 * (i + 1);
            y[i] = i + 1;
        }
        SimpleDataset data(x, y, yerr);
        
        LinearFitter fitter(data);
        REQUIRE(fitter.dof() == 8);
    }
}

