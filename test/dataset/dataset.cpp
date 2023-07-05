#include <catch2/catch_test_macros.hpp>

#include <dataset/Dataset.h>

#pragma once

#include <math/Matrix.h>
#include <dataset/PointSet.h>
#include <utility/Limit.h>
#include <io/ExistingFile.h>

TEST_CASE("Dataset::Dataset") {
    SECTION("default") {
        Dataset dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 0);
    }

    SECTION("Matrix&&") {
        Matrix<double> m(2, 2);
        Dataset dataset(std::move(m));
        CHECK(dataset.size() == 4);
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
        std::vector<std::string> col_names = {"a", "b", "c"};
        Dataset dataset(cols, col_names);
        CHECK(dataset.size() == 6);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.get_col_names() == col_names);
        CHECK(dataset.col("a") == Vector<double>{1, 4});
        CHECK(dataset.col("b") == Vector<double>{2, 5});
        CHECK(dataset.col("c") == Vector<double>{3, 6});
    }

    SECTION("vector<vector<double>>&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        Dataset dataset(cols);
        CHECK(dataset.size() == 6);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.col(0) == Vector<double>{1, 4});
        CHECK(dataset.col(1) == Vector<double>{2, 5});
        CHECK(dataset.col(2) == Vector<double>{3, 6});
    }

    SECTION("unsigned int, unsigned int") {
        Dataset dataset(2, 3);
        CHECK(dataset.size() == 6);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 3);
    }

    SECTION("ExistingFile&") {
        io::ExistingFile file("test/files/2epe.dat");
        Dataset dataset(file);
        CHECK(dataset.size() == 200);
        CHECK(dataset.size_rows() == 100);
        CHECK(dataset.size_cols() == 3);
    }
}

TEST_CASE("Dataset::col") {
    CHECK(false);
}

TEST_CASE("Dataset::row") {
    CHECK(false);
}

TEST_CASE("Dataset::size") {
    CHECK(false);
}

TEST_CASE("Dataset::empty") {
    CHECK(false);
}

TEST_CASE("Dataset::save") {
    CHECK(false);
}

TEST_CASE("Dataset::load") {
    CHECK(false);
}

TEST_CASE("Dataset::set_col_names") {
    CHECK(false);
}

TEST_CASE("Dataset::get_col_names") {
    CHECK(false);
}

TEST_CASE("Dataset::interpolate") {
    CHECK(false);
}

TEST_CASE("Dataset::rolling_average") {
    CHECK(false);
}

TEST_CASE("Dataset::append") {
    CHECK(false);
}

TEST_CASE("Dataset::limit_x") {
    CHECK(false);
}

TEST_CASE("Dataset::limit_y") {
    CHECK(false);
}

TEST_CASE("Dataset::x") {
    CHECK(false);
}

TEST_CASE("Dataset::y") {
    CHECK(false);
}