#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset2D.h>

using namespace ausaxs;

TEST_CASE("Dataset2D: stats") {
    Dataset2D data1(
        std::vector<double>{1, 1, 1, 1}, 
        std::vector<double>{10, 3, 5, 6}, 
        std::vector<double>{1, 2, 3, 4}
    );

    Dataset2D data2(
        std::vector<double>{1, 1, 1, 1, 1, 1}, 
        std::vector<double>{12, 14, 15, 15, 14, 17}, 
        std::vector<double>{1, 0.5, 0.1, 2, 0.9, 3}
    );

    Dataset2D data3(
        std::vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        std::vector<double>{54, 66, 78, 80, 82, 84, 84, 90, 93}, 
        std::vector<double>{1, 0.2, 2, 0.5, 0.9, 4, 1, 0.1, 0.4}
    );

    SECTION("weighted mean error") {
        CHECK_THAT(data1.weighted_mean_error(), Catch::Matchers::WithinAbs(0.838116, 1e-6));
        CHECK_THAT(data2.weighted_mean_error(), Catch::Matchers::WithinAbs(0.0968561, 1e-6));
        CHECK_THAT(data3.weighted_mean_error(), Catch::Matchers::WithinAbs(0.0848808, 1e-6));
    }

    SECTION("weighted mean") {
        CHECK_THAT(data1.weighted_mean(), Catch::Matchers::WithinAbs(8.204903, 1e-3));
        CHECK_THAT(data2.weighted_mean(), Catch::Matchers::WithinAbs(14.924834, 1e-3));
        CHECK_THAT(data3.weighted_mean(), Catch::Matchers::WithinAbs(85.125966, 1e-3));       
    }

    SECTION("mean") {
        CHECK_THAT(data1.mean(), Catch::Matchers::WithinAbs(6, 1e-3));
        CHECK_THAT(data2.mean(), Catch::Matchers::WithinAbs(14.5, 1e-3));
        CHECK_THAT(data3.mean(), Catch::Matchers::WithinAbs(79, 1e-3));
    }

    SECTION("std") {
        CHECK_THAT(data1.std(), Catch::Matchers::WithinAbs(2.943920, 1e-3));
        CHECK_THAT(data2.std(), Catch::Matchers::WithinAbs(1.643167, 1e-3));
        CHECK_THAT(data3.std(), Catch::Matchers::WithinAbs(12.103718, 1e-3));
    }
}
