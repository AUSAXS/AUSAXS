#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset.h>
#include <math/Matrix.h>
#include <dataset/PointSet.h>
#include <utility/Limit.h>
#include <io/ExistingFile.h>

struct fixture {
    std::vector<double> x = {   1,   2,   3,   4,   5,   6,   7,   8,   9};
    std::vector<double> y = {  -6,  -4,  -1,   2,   1,   3,   6,   7,   9};
    Dataset dataset = Dataset({x, y});
};

TEST_CASE("Dataset::Dataset") {
    SECTION("default") {
        Dataset dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 0);
    }

    SECTION("Dataset&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        std::vector<std::string> col_names = {"a", "b"};
        Dataset dataset(cols, col_names);
        Dataset dataset2(dataset);
        CHECK(dataset == dataset2);
    }

    SECTION("Dataset&&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        std::vector<std::string> col_names = {"a", "b"};
        Dataset dataset(cols, col_names);
        Dataset dataset2(std::move(dataset));
        REQUIRE(dataset2.size() == 3);
        REQUIRE(dataset2.size_rows() == 3);
        REQUIRE(dataset2.size_cols() == 2);
        CHECK(dataset2.get_col_names() == col_names);
        CHECK(dataset2.col("a") == Vector<double>{1, 2, 3});
        CHECK(dataset2.col("b") == Vector<double>{4, 5, 6});
    }

    SECTION("Matrix&&") {
        Matrix<double> m(2, 2);
        Dataset dataset(std::move(m));
        CHECK(dataset.size() == 2);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 2);
    }

    SECTION("vector<string>&") {
        std::vector<std::string> col_names = {"a", "b", "c"};
        Dataset dataset(col_names);
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.get_col_names() == col_names);
    }

    SECTION("vector<vector<double>>&, vector<string>&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        std::vector<std::string> col_names = {"a", "b"};
        Dataset dataset(cols, col_names);
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 2);
        CHECK(dataset.get_col_names() == col_names);
        CHECK(dataset.col("a") == Vector<double>{1, 2, 3});
        CHECK(dataset.col("b") == Vector<double>{4, 5, 6});
    }

    SECTION("vector<vector<double>>&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        Dataset dataset(cols);
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 2);
        CHECK(dataset.col(0) == Vector<double>{1, 2, 3});
        CHECK(dataset.col(1) == Vector<double>{4, 5, 6});
    }

    SECTION("unsigned int, unsigned int") {
        Dataset dataset(2, 3);
        CHECK(dataset.size() == 2);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 3);
    }

    SECTION("ExistingFile&") {
        io::ExistingFile file("test/files/2epe.dat");
        Dataset dataset(file);
        CHECK(dataset.size() == 104);
        CHECK(dataset.size_rows() == 104);
        CHECK(dataset.size_cols() == 4);
    }
}

TEST_CASE_METHOD(fixture, "Dataset::col") {
    CHECK(dataset.col(0) == Vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9});
    CHECK(dataset.col(1) == Vector<double>{-6, -4, -1, 2, 1, 3, 6, 7, 9});
}

TEST_CASE_METHOD(fixture, "Dataset::row") {
    CHECK(dataset.row(0) == Vector<double>{1, -6});
    CHECK(dataset.row(1) == Vector<double>{2, -4});
    CHECK(dataset.row(2) == Vector<double>{3, -1});
    CHECK(dataset.row(3) == Vector<double>{4, 2});
    CHECK(dataset.row(4) == Vector<double>{5, 1});
    CHECK(dataset.row(5) == Vector<double>{6, 3});
    CHECK(dataset.row(6) == Vector<double>{7, 6});
    CHECK(dataset.row(7) == Vector<double>{8, 7});
    CHECK(dataset.row(8) == Vector<double>{9, 9});
}

TEST_CASE_METHOD(fixture, "Dataset::size") {
    CHECK(dataset.size() == 9);
}

TEST_CASE_METHOD(fixture, "Dataset::empty") {
    CHECK(!dataset.empty());

    Dataset empty_dataset;
    CHECK(empty_dataset.empty());
}

TEST_CASE("Dataset::save") {
    Dataset dataset({
        std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10}, 
        std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10}
    });

    std::string path = "temp/dataset/save.dat";
    dataset.save(path);
    Dataset loaded_dataset(path);
    CHECK(dataset == loaded_dataset);
}

TEST_CASE_METHOD(fixture, "Dataset::col_names") {
    dataset.set_col_names({"a", "b"});
    CHECK(dataset.get_col_names() == std::vector<std::string>{"a", "b"});
}

TEST_CASE("Dataset::interpolate") {    
    SECTION("simple") {
        Dataset data({
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        });

        data = data.interpolate(1);
        REQUIRE(data.size() == 18);
        CHECK_THAT(data.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(data.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

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
        for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
            x.push_back(xx);
            y.push_back(sin(xx));
        }

        Dataset data({x, y});
        data = data.interpolate(5);
        for (unsigned int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(sin(data.x(i)), 1e-3));
        }
    }

    SECTION("vector interpolation") {
        std::vector<double> x1, y1, x2;
        for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
            x1.push_back(xx);
            y1.push_back(sin(xx));
            x2.push_back(xx + 0.025);
        }

        Dataset data1({x1, y1});
        auto data2 = data1.interpolate(x2);
        for (unsigned int i = 0; i < data2.size(); i++) {
            CHECK_THAT(data2.y(i), Catch::Matchers::WithinAbs(sin(data2.x(i)), 1e-3));
        }
    }

    SECTION("single values") {
        std::vector<double> x1, y1;
        for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
            x1.push_back(xx);
            y1.push_back(sin(xx));
        }

        Dataset data1({x1, y1});
        for (unsigned int i = 0; i < data1.size(); i++) {
            CHECK_THAT(data1.interpolate_y(data1.x(i)+0.025), Catch::Matchers::WithinAbs(sin(data1.x(i)+0.025), 1e-3));
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

TEST_CASE_METHOD(fixture, "Dataset::append") {
    auto d2 = dataset;
    dataset.append(d2);
    REQUIRE(dataset.size() == 18);
    for (unsigned int i = 0; i < 9; ++i) {
        CHECK(dataset.x(i) == dataset.x(i+9));
        CHECK(dataset.y(i) == dataset.y(i+9));
    }

    d2 = dataset;
    dataset.append(d2);
    REQUIRE(dataset.size() == 36);
    for (unsigned int i = 0; i < 18; ++i) {
        CHECK(dataset.x(i) == dataset.x(i+18));
        CHECK(dataset.y(i) == dataset.y(i+18));
    }
}

TEST_CASE_METHOD(fixture, "Dataset::limit_x") {
    SECTION("simple") {
        Limit limit(0.5, 5);
        dataset.limit_x(limit);
        for (unsigned int i = 0; i < dataset.size(); i++) {
            CHECK(limit.min <= dataset.x(i));
            CHECK(dataset.x(i) <= limit.max);
        }
    }

    SECTION("real data") {
        Dataset data("test/files/2epe.dat");

        unsigned int start = 0;
        while (data.x(start) < 0.01) {            
            start++;
        }

        unsigned int end = data.size()-1;
        while (0.3 < data.x(end)) {
            end--;
        }

        auto data_limited = data;
        data_limited.limit_x(0.01, 0.3);
        REQUIRE(data_limited.size() == end-start+1);
        for (unsigned int i = 0; i < data_limited.size(); i++) {
            CHECK(data_limited.x(i) == data.x(i+start));
            CHECK(data_limited.y(i) == data.y(i+start));
        }
    }
}

TEST_CASE_METHOD(fixture, "Dataset::limit_y") {
    Limit limit(0.5, 5);
    dataset.limit_y(limit);
    for (unsigned int i = 0; i < dataset.size(); i++) {
        CHECK(limit.min <= dataset.y(i));
        CHECK(dataset.y(i) <= limit.max);
    }
}

TEST_CASE_METHOD(fixture, "Dataset::x") {
    CHECK(dataset.x() == dataset.col(0));
}

TEST_CASE_METHOD(fixture, "Dataset::y") {
    CHECK(dataset.y() == dataset.col(1));
}

TEST_CASE("Dataset::sort_x") {
    std::vector<double> x = {1, 0, 4, 3, 2};
    std::vector<double> y = {10, 0, 40, 30, 20};
    Dataset data({x, y});

    data.sort_x();
    CHECK(data.x() == std::vector<double>({0, 1, 2, 3, 4}));
    CHECK(data.y() == std::vector<double>({0, 10, 20, 30, 40}));
}