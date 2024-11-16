#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset2D.h>

using namespace ausaxs;

TEST_CASE("Dataset2D_scaling_methods") {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {10, 20, 30, 40, 50};
    std::vector<double> yerr = {1, 2, 3, 4, 5};
    std::vector<double> xerr = {0.1, 0.2, 0.3, 0.4, 0.5};
    Dataset2D data(x, y, xerr, yerr);

    SECTION("scale y") {
        data.scale_y(2);
        CHECK(data.y() == std::vector<double>({20, 40, 60, 80, 100}));

        data.scale_y(0.2);
        CHECK(data.y() == std::vector<double>({4, 8, 12, 16, 20}));

        data.scale_y(-0.5);
        CHECK(data.y() == std::vector<double>({-2, -4, -6, -8, -10}));
    }

    SECTION("scale errors") {
        data.scale_errors(2);
        CHECK(data.yerr() == std::vector<double>({2, 4, 6, 8, 10}));
        CHECK(data.xerr() == std::vector<double>({0.2, 0.4, 0.6, 0.8, 1.0}));

        data.scale_errors(0.5);
        CHECK(data.yerr() == std::vector<double>({1, 2, 3, 4, 5}));
        CHECK(data.xerr() == std::vector<double>({0.1, 0.2, 0.3, 0.4, 0.5}));

        data.scale_errors(-0.5);
        CHECK(data.yerr() == std::vector<double>({-0.5, -1, -1.5, -2, -2.5}));
        CHECK(data.xerr() == std::vector<double>({-0.05, -0.1, -0.15, -0.2, -0.25}));
    }

    SECTION("normalize") {
        data.normalize(1);
        CHECK(data.y() == std::vector<double>({1, 2, 3, 4, 5}));

        data.normalize(-5);
        CHECK(data.y() == std::vector<double>({-5, -10, -15, -20, -25}));

        data.normalize(10);
        CHECK(data.y() == std::vector<double>({10, 20, 30, 40, 50}));
    }
}

TEST_CASE("stats") {
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

    SECTION("weighted_mean_error") {
        CHECK_THAT(data1.weighted_mean_error(), Catch::Matchers::WithinAbs(0.838116, 1e-6));
        CHECK_THAT(data2.weighted_mean_error(), Catch::Matchers::WithinAbs(0.0968561, 1e-6));
        CHECK_THAT(data3.weighted_mean_error(), Catch::Matchers::WithinAbs(0.0848808, 1e-6));
    }

    SECTION("weighted_mean") {
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

TEST_CASE("basics") {
    std::vector<double> xd = {   1,   2,   3,   4,   5,   6,   7,   8,   9};
    std::vector<double> yd = {  -6,  -4,  -1,   2,   1,   3,   6,   7,   9};
    std::vector<double> xed = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
    std::vector<double> yed = { -7,  -5,  -2,   1,   0,   2,   5,   6,   8};
    Dataset2D data(xd, yd, xed, yed);

    SECTION("columns") {
        auto x = data.x();
        auto y = data.y();
        auto yerr = data.yerr();
        auto xerr = data.xerr();

        REQUIRE(x.size() == 9);
        REQUIRE(y.size() == 9);
        REQUIRE(yerr.size() == 9);
        REQUIRE(xerr.size() == 9);
        CHECK(x == Vector(xd));
        CHECK(y == Vector(yd));
        CHECK(yerr == Vector(yed));
        CHECK(xerr == Vector(xed));

        // columns are mutable
        std::transform(x.begin(), x.end(), yerr.begin(), [](double x) {return x;});
        std::transform(y.begin(), y.end(), yerr.begin(), y.begin(), std::multiplies<>());
        CHECK(y == Vector{-6, -8, -3, 8, 5, 18, 42, 56, 81});
    }

    SECTION("column by name") {
        std::vector<std::string> names = {"x", "y", "yerr", "xerr"};
        data.set_col_names(names);
        auto x = data.col("x");
        auto y = data.col("y");
        auto yerr = data.col("yerr");
        auto xerr = data.col("xerr");

        REQUIRE(data.get_col_names() == names);
        CHECK(data.get_col_names(0) == "x");
        CHECK(data.get_col_names(1) == "y");
        CHECK(data.get_col_names(2) == "yerr");
        CHECK(data.get_col_names(3) == "xerr");

        CHECK(x == Vector(xd));
        CHECK(y == Vector(yd));
        CHECK(yerr == Vector(yed));
        CHECK(xerr == Vector(xed));
    }

    SECTION("indexing") {
        CHECK(data.x(0) == 1);
        CHECK(data.x(1) == 2);
        CHECK(data.x(2) == 3);

        CHECK(data.y(0) == -6);
        CHECK(data.y(1) == -4);
        CHECK(data.y(2) == -1);

        CHECK(data.xerr(0) == 0.5);
        CHECK(data.xerr(1) == 1.5);
        CHECK(data.xerr(2) == 2.5);

        CHECK(data.yerr(0) == -7);
        CHECK(data.yerr(1) == -5);
        CHECK(data.yerr(2) == -2);

        const auto& cdata = data;
        CHECK(cdata.x() == data.x());
        CHECK(cdata.x(0) == 1);
        CHECK(cdata.x(1) == 2);
        CHECK(cdata.x(2) == 3);

        CHECK(cdata.y() == data.y());
        CHECK(cdata.y(0) == -6);
        CHECK(cdata.y(1) == -4);
        CHECK(cdata.y(2) == -1);

        CHECK(cdata.xerr() == data.xerr());
        CHECK(cdata.xerr(0) == 0.5);
        CHECK(cdata.xerr(1) == 1.5);
        CHECK(cdata.xerr(2) == 2.5);

        CHECK(cdata.yerr() == data.yerr());
        CHECK(cdata.yerr(0) == -7);
        CHECK(cdata.yerr(1) == -5);
        CHECK(cdata.yerr(2) == -2);
    }

    SECTION("constructors") {
        Dataset2D data2(data);
        CHECK(data2.x() == data.x());
        CHECK(data2.y() == data.y());
        CHECK(data2.xerr() == data.xerr());
        CHECK(data2.yerr() == data.yerr());

        Dataset2D data3(std::move(data));
        CHECK(data3.x() == data2.x());
        CHECK(data3.y() == data2.y());
        CHECK(data3.xerr() == data2.xerr());
        CHECK(data3.yerr() == data2.yerr());

        Dataset data4(std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
        CHECK(data4.x() == Vector{1, 2, 3});
        CHECK(data4.y() == Vector{4, 5, 6});
        CHECK(data4.col(2) == Vector{7, 8, 9});
        SimpleDataset data5(data4);
        CHECK(data5.x() == Vector{1, 2, 3});
        CHECK(data5.y() == Vector{4, 5, 6});
        CHECK(data5.yerr() == Vector{7, 8, 9});
        CHECK(data5.N == 3);
        CHECK(data5.M == 3);

        data4 = Dataset(std::vector<std::vector<double>>{{1, 3, 5, 7, 9}, {2, 4, 6, 8, 10}});
        CHECK(data4.x() == Vector{1, 3, 5, 7, 9});
        CHECK(data4.y() == Vector{2, 4, 6, 8, 10});
        SimpleDataset data6(data4);
        CHECK(data6.x() == Vector{1, 3, 5, 7, 9});
        CHECK(data6.y() == Vector{2, 4, 6, 8, 10});
        CHECK(data6.N == 5);
        CHECK(data6.M == 3);
    }
}


TEST_CASE("Dataset2D::pushback") {
    SECTION("various") {
        Dataset2D data;
        data.push_back(1, 2, 3, 4);
        data.push_back(5, 6, 7, 8);
        data.push_back(9, 10, 11, 12);
        data.push_back(Point2D(13, 14, 15, 16));        

        CHECK(data.x() == Vector{1, 5, 9, 13});
        CHECK(data.y() == Vector{2, 6, 10, 14});
        CHECK(data.xerr() == Vector{3, 7, 11, 15});
        CHECK(data.yerr() == Vector{4, 8, 12, 16});
    }

    SECTION("Point2D") {
        Point2D p0(0, 4);
        Point2D p1(1, 2);
        Point2D p2(3, 4);
        Point2D p3(5, 6);
        Point2D p4(7, 8);
        Point2D p5(9, 10);

        Dataset2D data;
        data.push_back(p0);
        data.push_back(p1);
        data.push_back(p2);
        data.push_back(p3);
        data.push_back(p4);
        data.push_back(p5);

        CHECK(data.size() == 6);
        CHECK(data.x() == std::vector<double>({0, 1, 3, 5, 7, 9}));
        CHECK(data.y() == std::vector<double>({4, 2, 4, 6, 8, 10}));

        CHECK(data.get_point(0) == p0);
        CHECK(data.get_point(1) == p1);
        CHECK(data.get_point(2) == p2);

        CHECK(data.find_minimum() == p1);
    }
}
