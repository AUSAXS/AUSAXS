#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset2D.h>

using namespace ausaxs;

TEST_CASE("Dataset2D::Dataset2D") {
    SECTION("default constructor") {
        Dataset2D dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 4);
    }

    SECTION("unsigned int") {
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

    SECTION("vector<double>, vector<double>, string, string") {
        Dataset2D dataset({1, 2, 3}, {4, 5, 6}, "xlabel", "ylabel");
        CHECK(dataset.size() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.get_col_names()[0] == "xlabel");
        CHECK(dataset.get_col_names()[1] == "ylabel");
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
    std::vector<double> xv = {1, 2, 3};
    std::vector<double> yv = {4, 5, 6};
    std::vector<double> xev = {7, 8, 9};
    std::vector<double> yev = {10, 11, 12};
    Dataset2D dataset(xv, yv, xev, yev);

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
    std::vector<double> xv = {1, 2, 3};
    std::vector<double> yv = {10, 20, 30};
    std::vector<double> xev = {1, 2, 3};
    std::vector<double> yev = {4, 5, 6};
    Dataset2D dataset(xv, yv, xev, yev);
    
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
    CHECK(x == Vector(xd));
    CHECK(y == Vector(yd));
    CHECK(yerr == Vector(yed));
    CHECK(xerr == Vector(xed));
}

TEST_CASE("Dataset2D::indexing") {
    std::vector<double> xd = {1, 2, 3};
    std::vector<double> yd = {4, 5, 6};
    std::vector<double> xed = {0.1, 0.2, 0.3};
    std::vector<double> yed = {7, 8, 9};
    Dataset2D data(xd, yd, xed, yed);

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
