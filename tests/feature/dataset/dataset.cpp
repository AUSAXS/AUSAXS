#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset.h>
#include <io/ExistingFile.h>
#include <settings/GeneralSettings.h>

#include <support/temp_file.h>

#include <numbers>

using namespace ausaxs;

struct fixture {
    std::vector<double> x = {   1,   2,   3,   4,   5,   6,   7,   8,   9};
    std::vector<double> y = {  -6,  -4,  -1,   2,   1,   3,   6,   7,   9};
    Dataset dataset = Dataset({x, y});
};

TEST_CASE("Dataset::Dataset") {
    settings::general::verbose = false;
    SECTION("ExistingFile&") {
        io::ExistingFile file("tests/files/2epe.dat");
        Dataset dataset(file);
        CHECK(dataset.size() == 104);
        CHECK(dataset.size_rows() == 104);
        CHECK(dataset.size_cols() == 4);
    }
}

TEST_CASE("Dataset::save") {
    settings::general::verbose = false;

    Dataset dataset({
        std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10}, 
        std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10}
    });

    test::TempFile path(".dat");
    dataset.save(path);
    Dataset loaded_dataset(path);
    CHECK(dataset == loaded_dataset);
}



TEST_CASE("Dataset::interpolate") {    
    SECTION("simple") {
        Dataset data({
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        });

        data = data.interpolate(1);
        REQUIRE(data.size() == 19); // 10 originals + 1 inserted in each of the 9 gaps
        CHECK_THAT(data.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(data.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

        // the final point of the input must survive the interpolation
        CHECK_THAT(data.x(18), Catch::Matchers::WithinAbs(10, 1e-6));
        CHECK_THAT(data.y(18), Catch::Matchers::WithinAbs(10, 1e-6));

        CHECK_THAT(data.x(17), Catch::Matchers::WithinAbs(9.5, 1e-6));
        CHECK_THAT(data.y(17), Catch::Matchers::WithinAbs(9.5, 1e-6));

        CHECK_THAT(data.x(1), Catch::Matchers::WithinAbs(1.5, 1e-6));
        CHECK_THAT(data.y(1), Catch::Matchers::WithinAbs(1.5, 1e-6));

        CHECK_THAT(data.x(2), Catch::Matchers::WithinAbs(2, 1e-6));
        CHECK_THAT(data.y(2), Catch::Matchers::WithinAbs(2, 1e-6));

        CHECK_THAT(data.x(3), Catch::Matchers::WithinAbs(2.5, 1e-6));
        CHECK_THAT(data.y(3), Catch::Matchers::WithinAbs(2.5, 1e-6));

        CHECK_THAT(data.x(4), Catch::Matchers::WithinAbs(3, 1e-6));
        CHECK_THAT(data.y(4), Catch::Matchers::WithinAbs(3, 1e-6));

        CHECK_THAT(data.x(5), Catch::Matchers::WithinAbs(3.5, 1e-6));
        CHECK_THAT(data.y(5), Catch::Matchers::WithinAbs(3.5, 1e-6));
    }

    SECTION("sine") {
        std::vector<double> x, y;
        for (double xx = 0; xx < 2*std::numbers::pi; xx += 0.05) {
            x.push_back(xx);
            y.push_back(std::sin(xx));
        }

        Dataset data({x, y});
        data = data.interpolate(5);
        for (int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(std::sin(data.x(i)), 1e-3));
        }
    }

    SECTION("vector interpolation") {
        std::vector<double> x1, y1, x2;
        for (double xx = 0; xx < 2*std::numbers::pi; xx += 0.05) {
            x1.push_back(xx);
            y1.push_back(std::sin(xx));
            x2.push_back(xx + 0.025);
        }

        Dataset data1({x1, y1});
        auto data2 = data1.interpolate(x2);
        for (int i = 0; i < data2.size(); i++) {
            CHECK_THAT(data2.y(i), Catch::Matchers::WithinAbs(std::sin(data2.x(i)), 1e-3));
        }
    }

    SECTION("single values") {
        std::vector<double> x1, y1;
        for (double xx = 0; xx < 2*std::numbers::pi; xx += 0.05) {
            x1.push_back(xx);
            y1.push_back(std::sin(xx));
        }

        Dataset data1({x1, y1});
        for (int i = 0; i < data1.size(); i++) {
            CHECK_THAT(data1.interpolate_x(data1.x(i)+0.025, 1), Catch::Matchers::WithinAbs(std::sin(data1.x(i)+0.025), 1e-3));
        }
    }

    SECTION("multiple columns") {
        std::vector<double> x1, y1, x2, y2;
        for (double x = 0; x < 2*std::numbers::pi; x += 0.05) {
            x1.push_back(x);
            x2.push_back(x + 0.025);
            y1.push_back(std::sin(x));
            y2.push_back(std::cos(x));
        }

        Dataset data({x1, y1, y2});
        auto data3 = data.interpolate(x2);
        REQUIRE(data3.x() == x2);
        for (int i = 0; i < data3.size(); i++) {
            CHECK_THAT(data3.col(1)[i], Catch::Matchers::WithinAbs(std::sin(data3.x(i)), 1e-3));
            CHECK_THAT(data3.col(2)[i], Catch::Matchers::WithinAbs(std::cos(data3.x(i)), 1e-3));
            CHECK_THAT(data.interpolate_x(data3.x(i), 1), Catch::Matchers::WithinAbs(std::sin(data3.x(i)), 1e-3));
            CHECK_THAT(data.interpolate_x(data3.x(i), 2), Catch::Matchers::WithinAbs(std::cos(data3.x(i)), 1e-3));
        }
    }
}

TEST_CASE("Dataset::rolling_average") {
    Dataset data({
        std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
        std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
        std::vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
    });

    SECTION("half_moving_average") {
        SECTION("3") {
            Dataset res = data.rolling_average(3);
            REQUIRE(res.size() == 10);
            CHECK_THAT(res.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
            CHECK_THAT(res.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

            CHECK_THAT(res.x(1), Catch::Matchers::WithinAbs(2, 1e-6));
            CHECK_THAT(res.y(1), Catch::Matchers::WithinAbs((1./2 + 2 + 3./2)/2, 1e-6));

            CHECK_THAT(res.x(2), Catch::Matchers::WithinAbs(3, 1e-6));
            CHECK_THAT(res.y(2), Catch::Matchers::WithinAbs((2./2 + 3 + 4./2)/2, 1e-6));

            CHECK_THAT(res.x(3), Catch::Matchers::WithinAbs(4, 1e-6));
            CHECK_THAT(res.y(3), Catch::Matchers::WithinAbs((3./2 + 4 + 5./2)/2, 1e-6));

            CHECK_THAT(res.x(4), Catch::Matchers::WithinAbs(5, 1e-6));
            CHECK_THAT(res.y(4), Catch::Matchers::WithinAbs((4./2 + 5 + 6./2)/2, 1e-6));

            CHECK_THAT(res.x(5), Catch::Matchers::WithinAbs(6, 1e-6));
            CHECK_THAT(res.y(5), Catch::Matchers::WithinAbs((5./2 + 6 + 7./2)/2, 1e-6));

            CHECK_THAT(res.x(6), Catch::Matchers::WithinAbs(7, 1e-6));
            CHECK_THAT(res.y(6), Catch::Matchers::WithinAbs((6./2 + 7 + 8./2)/2, 1e-6));

            CHECK_THAT(res.x(7), Catch::Matchers::WithinAbs(8, 1e-6));
            CHECK_THAT(res.y(7), Catch::Matchers::WithinAbs((7./2 + 8 + 9./2)/2, 1e-6));

            CHECK_THAT(res.x(8), Catch::Matchers::WithinAbs(9, 1e-6));
            CHECK_THAT(res.y(8), Catch::Matchers::WithinAbs((8./2 + 9 + 10./2)/2, 1e-6));

            CHECK_THAT(res.x(9), Catch::Matchers::WithinAbs(10, 1e-6));
            CHECK_THAT(res.y(9), Catch::Matchers::WithinAbs(10, 1e-6));
        }

        SECTION("5") {
            Dataset res = data.rolling_average(5);
            REQUIRE(res.size() == 10);
            CHECK_THAT(res.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
            CHECK_THAT(res.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

            CHECK_THAT(res.x(1), Catch::Matchers::WithinAbs(2, 1e-6));
            CHECK_THAT(res.y(1), Catch::Matchers::WithinAbs((1./2 + 2 + 3./2)/2, 1e-6));

            CHECK_THAT(res.x(2), Catch::Matchers::WithinAbs(3, 1e-6));
            CHECK_THAT(res.y(2), Catch::Matchers::WithinAbs((1./4 + 2./2 + 3 + 4./2 + 5./4)/2.5, 1e-6));

            CHECK_THAT(res.x(3), Catch::Matchers::WithinAbs(4, 1e-6));
            CHECK_THAT(res.y(3), Catch::Matchers::WithinAbs((2./4 + 3./2 + 4 + 5./2 + 6./4)/2.5, 1e-6));

            CHECK_THAT(res.x(4), Catch::Matchers::WithinAbs(5, 1e-6));
            CHECK_THAT(res.y(4), Catch::Matchers::WithinAbs((3./4 + 4./2 + 5 + 6./2 + 7./4)/2.5, 1e-6));

            CHECK_THAT(res.x(5), Catch::Matchers::WithinAbs(6, 1e-6));
            CHECK_THAT(res.y(5), Catch::Matchers::WithinAbs((4./4 + 5./2 + 6 + 7./2 + 8./4)/2.5, 1e-6));

            CHECK_THAT(res.x(6), Catch::Matchers::WithinAbs(7, 1e-6));
            CHECK_THAT(res.y(6), Catch::Matchers::WithinAbs((5./4 + 6./2 + 7 + 8./2 + 9./4)/2.5, 1e-6));

            CHECK_THAT(res.x(7), Catch::Matchers::WithinAbs(8, 1e-6));
            CHECK_THAT(res.y(7), Catch::Matchers::WithinAbs((6./4 + 7./2 + 8 + 9./2 + 10./4)/2.5, 1e-6));

            CHECK_THAT(res.x(8), Catch::Matchers::WithinAbs(9, 1e-6));
            CHECK_THAT(res.y(8), Catch::Matchers::WithinAbs((8./2 + 9 + 10./2)/2, 1e-6));

            CHECK_THAT(res.x(9), Catch::Matchers::WithinAbs(10, 1e-6));
            CHECK_THAT(res.y(9), Catch::Matchers::WithinAbs(10, 1e-6));
        }
    }
}

TEST_CASE_METHOD(fixture, "Dataset::limit") {
    settings::general::verbose = false;
    SECTION("real data") {
        Dataset data("tests/files/2epe.dat");

        int start = 0;
        while (data.x(start) < 0.01) {            
            start++;
        }

        int end = data.size()-1;
        while (0.3 < data.x(end)) {
            end--;
        }

        auto data_limited = data;
        data_limited.limit(0, 0.01, 0.3);
        REQUIRE(data_limited.size() == end-start+1);
        for (int i = 0; i < data_limited.size(); i++) {
            CHECK(data_limited.x(i) == data.x(i+start));
            CHECK(data_limited.y(i) == data.y(i+start));
        }
    }
}
