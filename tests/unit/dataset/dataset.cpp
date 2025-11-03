#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset.h>
#include <math/Matrix.h>
#include <utility/Limit.h>

using namespace ausaxs;

TEST_CASE("Dataset::Dataset") {
    SECTION("default constructor") {
        Dataset dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.empty());
    }

    SECTION("copy constructor") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        Dataset dataset(cols);
        Dataset dataset2(dataset);
        CHECK(dataset == dataset2);
        CHECK(dataset2.col(0) == Vector<double>{1, 2, 3});
        CHECK(dataset2.col(1) == Vector<double>{4, 5, 6});
    }

    SECTION("move constructor") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        Dataset dataset(cols);
        Dataset dataset2(std::move(dataset));
        REQUIRE(dataset2.size() == 3);
        REQUIRE(dataset2.size_rows() == 3);
        REQUIRE(dataset2.size_cols() == 2);
        CHECK(dataset2.col(0) == Vector<double>{1, 2, 3});
        CHECK(dataset2.col(1) == Vector<double>{4, 5, 6});
    }

    SECTION("Matrix&&") {
        Matrix<double> m(2, 2);
        Dataset dataset(std::move(m));
        CHECK(dataset.size() == 2);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 2);
    }

    SECTION("vector<vector<double>>&, vector<string>&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        Dataset dataset(cols);
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 2);
        CHECK(dataset.col(0) == Vector<double>{1, 2, 3});
        CHECK(dataset.col(1) == Vector<double>{4, 5, 6});
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
}

TEST_CASE("Dataset::col") {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {6, 7, 8, 9, 10};
    Dataset dataset({x, y});

    SECTION("by index") {
        CHECK(dataset.col(0) == Vector<double>{1, 2, 3, 4, 5});
        CHECK(dataset.col(1) == Vector<double>{6, 7, 8, 9, 10});
    }
}

TEST_CASE("Dataset::row") {
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    Dataset dataset({x, y});

    CHECK(dataset.row(0) == Vector<double>{1, 4});
    CHECK(dataset.row(1) == Vector<double>{2, 5});
    CHECK(dataset.row(2) == Vector<double>{3, 6});
}

TEST_CASE("Dataset::x") {
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    Dataset dataset({x, y});

    SECTION("column accessor") {
        auto _x = dataset.x();
        CHECK(_x == dataset.col(0));
    }

    SECTION("element accessor") {
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.x(1) == 2);
        CHECK(dataset.x(2) == 3);
    }
}

TEST_CASE("Dataset::y") {
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    Dataset dataset({x, y});

    SECTION("column accessor") {
        CHECK(dataset.y() == dataset.col(1));
    }

    SECTION("element accessor") {
        CHECK(dataset.y(0) == 4);
        CHECK(dataset.y(1) == 5);
        CHECK(dataset.y(2) == 6);
    }
}

TEST_CASE("Dataset::size") {
    SECTION("empty dataset") {
        Dataset dataset;
        CHECK(dataset.size() == 0);
    }

    SECTION("non-empty dataset") {
        Dataset dataset(5, 3);
        CHECK(dataset.size() == 5);
        CHECK(dataset.size_rows() == 5);
        CHECK(dataset.size_cols() == 3);
    }
}

TEST_CASE("Dataset::empty") {
    SECTION("empty dataset") {
        Dataset dataset;
        CHECK(dataset.empty());
    }

    SECTION("non-empty dataset") {
        Dataset dataset(5, 3);
        CHECK_FALSE(dataset.empty());
    }
}

TEST_CASE("Dataset::push_back") {
    Dataset dataset(0, 3);
    CHECK(dataset.size() == 0);

    dataset.push_back({1, 2, 3});
    CHECK(dataset.size() == 1);
    CHECK(dataset.row(0) == Vector<double>{1, 2, 3});

    dataset.push_back({4, 5, 6});
    CHECK(dataset.size() == 2);
    CHECK(dataset.row(1) == Vector<double>{4, 5, 6});
}

TEST_CASE("Dataset::operator==") {
    std::vector<std::vector<double>> cols1 = {{1, 2, 3}, {4, 5, 6}};
    std::vector<std::vector<double>> cols2 = {{1, 2, 3}, {4, 5, 7}};
    
    Dataset dataset1(cols1);
    Dataset dataset2(cols1);
    Dataset dataset3(cols2);

    CHECK(dataset1 == dataset2);
    CHECK_FALSE(dataset1 == dataset3);
}

TEST_CASE("Dataset::index") {
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    Dataset dataset({x, y});

    CHECK(dataset.index(0, 0) == 1);
    CHECK(dataset.index(0, 1) == 4);
    CHECK(dataset.index(1, 0) == 2);
    CHECK(dataset.index(1, 1) == 5);
    CHECK(dataset.index(2, 0) == 3);
    CHECK(dataset.index(2, 1) == 6);
}

TEST_CASE("Dataset::select_columns") {
    std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Dataset dataset(cols);

    SECTION("by index") {
        Dataset selected = dataset.select_columns({0, 2});
        CHECK(selected.size_cols() == 2);
        CHECK(selected.col(0) == Vector<double>{1, 2, 3});
        CHECK(selected.col(1) == Vector<double>{7, 8, 9});
    }
}

TEST_CASE("Dataset::limit_x") {
    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> y = {10, 20, 30, 40, 50, 60, 70, 80, 90};
    Dataset dataset({x, y});

    dataset.limit_x(3, 7);
    CHECK(dataset.size() == 5);
    CHECK(dataset.x(0) == 3);
    CHECK(dataset.x(4) == 7);
}

TEST_CASE("Dataset::limit_y") {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {10, 20, 30, 40, 50};
    Dataset dataset({x, y});

    dataset.limit_y(20, 40);
    CHECK(dataset.size() == 3);
    CHECK(dataset.y(0) == 20);
    CHECK(dataset.y(1) == 30);
    CHECK(dataset.y(2) == 40);
}

TEST_CASE("Dataset::append") {
    std::vector<double> x1 = {1, 2, 3};
    std::vector<double> y1 = {4, 5, 6};
    Dataset dataset1({x1, y1});

    std::vector<double> x2 = {7, 8};
    std::vector<double> y2 = {9, 10};
    Dataset dataset2({x2, y2});

    dataset1.append(dataset2);
    CHECK(dataset1.size() == 5);
    CHECK(dataset1.x(3) == 7);
    CHECK(dataset1.x(4) == 8);
    CHECK(dataset1.y(3) == 9);
    CHECK(dataset1.y(4) == 10);
}

TEST_CASE("Dataset::sort_x") {
    std::vector<double> x = {5, 2, 8, 1, 9};
    std::vector<double> y = {10, 20, 30, 40, 50};
    Dataset dataset({x, y});

    dataset.sort_x();
    CHECK(dataset.x(0) == 1);
    CHECK(dataset.x(1) == 2);
    CHECK(dataset.x(2) == 5);
    CHECK(dataset.x(3) == 8);
    CHECK(dataset.x(4) == 9);
    CHECK(dataset.y(0) == 40);
    CHECK(dataset.y(1) == 20);
    CHECK(dataset.y(2) == 10);
    CHECK(dataset.y(3) == 30);
    CHECK(dataset.y(4) == 50);
}
