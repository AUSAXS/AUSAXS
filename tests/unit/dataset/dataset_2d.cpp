#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset2D.h>

using namespace ausaxs;

TEST_CASE("Dataset2D::Dataset2D") {
    SECTION("default constructor") {
        Dataset2D dataset;
        CHECK(dataset.empty());
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 4);
    }

    SECTION("int") {
        Dataset2D dataset(10);
        CHECK(dataset.size() == 10);
        CHECK(dataset.size_rows() == 10);
        CHECK(dataset.size_cols() == 4);
    }

    SECTION("vector<double>, vector<double>") {
        Dataset2D dataset({1, 2, 3}, {4, 5, 6});
        CHECK(dataset.size() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{0, 0, 0});
        CHECK(dataset.xerr() == std::vector{0, 0, 0});
    }

    SECTION("vector<double>, vector<double>, vector<double>") {
        Dataset2D dataset({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
        CHECK(dataset.size() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{7, 8, 9});
        CHECK(dataset.xerr() == std::vector{0, 0, 0});
    }

    SECTION("vector<double>, vector<double>, vector<double>, vector<double>") {
        std::vector<double> xv = {1, 2, 3};
        std::vector<double> yv = {4, 5, 6};
        std::vector<double> xev = {7, 8, 9};
        std::vector<double> yev = {10, 11, 12};
        Dataset2D dataset(xv, yv, xev, yev);
        CHECK(dataset.size() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.xerr() == std::vector{7, 8, 9});
        CHECK(dataset.yerr() == std::vector{10, 11, 12});
    }

    SECTION("SimpleDataset") {
        SimpleDataset simple({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
        Dataset2D dataset(simple);
        CHECK(dataset.size() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{7, 8, 9});
    }
}

TEST_CASE("Dataset2D::xerr") {
    Dataset2D dataset({1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12});

    SECTION("column accessor") {
        CHECK(dataset.xerr() == dataset.col(3));
    }

    SECTION("element accessor") {
        CHECK(dataset.xerr(0) == 7);
        CHECK(dataset.xerr(1) == 8);
        CHECK(dataset.xerr(2) == 9);
    }
}

TEST_CASE("Dataset2D::push_back") {
    Dataset2D dataset;

    SECTION("x, y, xerr, yerr") {
        dataset.push_back(1, 2, 3, 4);
        CHECK(dataset.size() == 1);
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.y(0) == 2);
        CHECK(dataset.xerr(0) == 3);
        CHECK(dataset.yerr(0) == 4);

        dataset.push_back(5, 6, 7, 8);
        CHECK(dataset.size() == 2);
        CHECK(dataset.x(1) == 5);
        CHECK(dataset.y(1) == 6);
        CHECK(dataset.xerr(1) == 7);
        CHECK(dataset.yerr(1) == 8);
    }

    SECTION("x, y") {
        dataset.push_back(1, 2);
        CHECK(dataset.size() == 1);
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.y(0) == 2);
        CHECK(dataset.xerr(0) == 0);
        CHECK(dataset.yerr(0) == 0);
    }

    SECTION("Point2D") {
        Point2D point1(1, 2, 3, 4);
        dataset.push_back(point1);
        CHECK(dataset.size() == 1);
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.y(0) == 2);
        CHECK(dataset.xerr(0) == 3);
        CHECK(dataset.yerr(0) == 4);
    }
}

TEST_CASE("Dataset2D::scale_errors") {
    Dataset2D dataset({1, 2, 3}, {10, 20, 30}, {1, 2, 3}, {4, 5, 6});
    
    dataset.scale_errors(2);
    CHECK(dataset.y(0) == 10);
    CHECK(dataset.y(1) == 20);
    CHECK(dataset.y(2) == 30);
    CHECK(dataset.xerr(0) == 2);
    CHECK(dataset.xerr(1) == 4);
    CHECK(dataset.xerr(2) == 6);
    CHECK(dataset.yerr(0) == 8);
    CHECK(dataset.yerr(1) == 10);
    CHECK(dataset.yerr(2) == 12);
}

TEST_CASE("Dataset2D::columns") {
    std::vector<double> xd = {1, 2, 3, 4, 5};
    std::vector<double> yd = {10, 20, 30, 40, 50};
    std::vector<double> xed = {0.1, 0.2, 0.3, 0.4, 0.5};
    std::vector<double> yed = {1, 2, 3, 4, 5};
    Dataset2D data(xd, yd, xed, yed);

    auto x = data.x();
    auto y = data.y();
    auto yerr = data.yerr();
    auto xerr = data.xerr();

    REQUIRE(x.size() == 5);
    REQUIRE(y.size() == 5);
    REQUIRE(yerr.size() == 5);
    REQUIRE(xerr.size() == 5);
    CHECK(x == xd);
    CHECK(y == yd);
    CHECK(yerr == yed);
    CHECK(xerr == xed);
}

TEST_CASE("Dataset2D::indexing") {
    Dataset2D data({1, 2, 3}, {4, 5, 6}, {0.1, 0.2, 0.3}, {7, 8, 9});

    CHECK(data.x(0) == 1);
    CHECK(data.x(1) == 2);
    CHECK(data.x(2) == 3);

    CHECK(data.y(0) == 4);
    CHECK(data.y(1) == 5);
    CHECK(data.y(2) == 6);

    CHECK(data.xerr(0) == 0.1);
    CHECK(data.xerr(1) == 0.2);
    CHECK(data.xerr(2) == 0.3);

    CHECK(data.yerr(0) == 7);
    CHECK(data.yerr(1) == 8);
    CHECK(data.yerr(2) == 9);

    SECTION("const accessors") {
        const auto& cdata = data;
        CHECK(cdata.x(0) == 1);
        CHECK(cdata.y(0) == 4);
        CHECK(cdata.xerr(0) == 0.1);
        CHECK(cdata.yerr(0) == 7);
    }
}

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
        CHECK_THAT(data1.mean(1), Catch::Matchers::WithinAbs(6, 1e-3));
        CHECK_THAT(data2.mean(1), Catch::Matchers::WithinAbs(14.5, 1e-3));
        CHECK_THAT(data3.mean(1), Catch::Matchers::WithinAbs(79, 1e-3));
    }

    SECTION("std") {
        CHECK_THAT(data1.std(1), Catch::Matchers::WithinAbs(2.943920, 1e-3));
        CHECK_THAT(data2.std(1), Catch::Matchers::WithinAbs(1.643167, 1e-3));
        CHECK_THAT(data3.std(1), Catch::Matchers::WithinAbs(12.103718, 1e-3));
    }
}
