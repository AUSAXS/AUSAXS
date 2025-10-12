#include <catch2/catch_test_macros.hpp>

#include <math/MovingAverager.h>
#include <dataset/SimpleDataset.h>

using namespace ausaxs;

TEST_CASE("moving_average") {
    SECTION("simple") {
        SimpleDataset data(
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            std::vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
        );
        SECTION("unweighted") {
            SECTION("3") {
                auto res = MovingAverage::average(data.y(), 3);
                REQUIRE(res.size() == 10);
                CHECK(res[0] == 1);
                CHECK(res[1] == 2);
                CHECK(res[2] == 3);
                CHECK(res[3] == 4);
                CHECK(res[4] == 5);
                CHECK(res[5] == 6);
                CHECK(res[6] == 7);
                CHECK(res[7] == 8);
                CHECK(res[8] == 9);
                CHECK(res[9] == 10);
            }

            SECTION("5") {
                auto res = MovingAverage::average(data.y(), 5);
                REQUIRE(res.size() == 10);
                CHECK(res[0] == 1);
                CHECK(res[1] == 2);
                CHECK(res[2] == 3);
                CHECK(res[3] == 4);
                CHECK(res[4] == 5);
                CHECK(res[5] == 6);
                CHECK(res[6] == 7);
                CHECK(res[7] == 8);
                CHECK(res[8] == 9);
                CHECK(res[9] == 10);
            }
        }
    }
}