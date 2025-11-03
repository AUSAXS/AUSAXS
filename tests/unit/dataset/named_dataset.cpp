// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <dataset/NamedDataset.h>

using namespace ausaxs;

TEST_CASE("NamedDataset::NamedDataset") {
    SECTION("default constructor") {
        NamedDataset dataset;
        CHECK(dataset.empty());
        CHECK(dataset.names.empty());
    }

    SECTION("from Dataset with default names") {
        Dataset d(10, 3);
        NamedDataset named(std::move(d));
        CHECK(named.size_rows() == 10);
        CHECK(named.size_cols() == 3);
        CHECK(named.names.size() == 3);
        CHECK(named.names[0] == "col_0");
        CHECK(named.names[1] == "col_1");
        CHECK(named.names[2] == "col_2");
    }

    SECTION("from Dataset with explicit names") {
        Dataset d(10, 3);
        std::vector<std::string> names = {"x", "y", "z"};
        NamedDataset named(std::move(d), names);
        CHECK(named.size_rows() == 10);
        CHECK(named.size_cols() == 3);
        CHECK(named.names == names);
    }

    SECTION("from SimpleDataset") {
        SimpleDataset sd({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {0.1, 0.2, 0.3});
        NamedSimpleDataset named(std::move(sd), {"q", "I", "I_err"});
        CHECK(named.size() == 3);
        CHECK(named.names == std::vector<std::string>{"q", "I", "I_err"});
    }
}

TEST_CASE("NamedDataset::set_col_names") {
    Dataset d(10, 3);
    NamedDataset named(std::move(d));

    SECTION("set all names") {
        std::vector<std::string> names = {"first", "second", "third"};
        named.set_col_names(names);
        CHECK(named.names == names);
    }

    SECTION("set individual name") {
        named.set_col_names(0, "first");
        named.set_col_names(1, "second");
        named.set_col_names(2, "third");
        CHECK(named.names[0] == "first");
        CHECK(named.names[1] == "second");
        CHECK(named.names[2] == "third");
    }

    SECTION("mismatched size throws") {
        std::vector<std::string> names = {"first", "second"};
        CHECK_THROWS(named.set_col_names(names));
    }
}

TEST_CASE("NamedDataset::get_col_names") {
    Dataset d(10, 3);
    std::vector<std::string> names = {"x", "y", "z"};
    NamedDataset named(std::move(d), names);

    SECTION("get all names") {
        auto result = named.get_col_names();
        CHECK(result == names);
    }

    SECTION("get individual name") {
        CHECK(named.get_col_names(0) == "x");
        CHECK(named.get_col_names(1) == "y");
        CHECK(named.get_col_names(2) == "z");
    }
}

TEST_CASE("NamedDataset::is_named") {
    SECTION("default names") {
        Dataset d(10, 3);
        NamedDataset named(std::move(d));
        CHECK_FALSE(named.is_named());
    }

    SECTION("custom names") {
        Dataset d(10, 3);
        NamedDataset named(std::move(d), {"x", "y", "z"});
        CHECK(named.is_named());
    }

    SECTION("partially custom names") {
        Dataset d(10, 3);
        NamedDataset named(std::move(d));
        named.set_col_names(1, "custom");
        CHECK(named.is_named());
    }
}

TEST_CASE("NamedDataset::col") {
    std::vector<double> x_data = {1.0, 2.0, 3.0};
    std::vector<double> y_data = {4.0, 5.0, 6.0};
    std::vector<double> z_data = {7.0, 8.0, 9.0};
    Dataset d({x_data, y_data, z_data});
    NamedDataset named(std::move(d), {"x", "y", "z"});

    SECTION("access by name") {
        auto col_x = named.col("x");
        CHECK(col_x[0] == 1.0);
        CHECK(col_x[1] == 2.0);
        CHECK(col_x[2] == 3.0);

        auto col_y = named.col("y");
        CHECK(col_y[0] == 4.0);
        CHECK(col_y[1] == 5.0);
        CHECK(col_y[2] == 6.0);
    }

    SECTION("invalid name throws") {
        CHECK_THROWS(named.col("invalid"));
    }
}

TEST_CASE("NamedDataset::select_columns") {
    std::vector<double> x_data = {1.0, 2.0, 3.0};
    std::vector<double> y_data = {4.0, 5.0, 6.0};
    std::vector<double> z_data = {7.0, 8.0, 9.0};
    Dataset d({x_data, y_data, z_data});
    NamedDataset named(std::move(d), {"x", "y", "z"});

    SECTION("select by name") {
        auto selected = named.select_columns({"x", "z"});
        CHECK(selected.size_cols() == 2);
        CHECK(selected.names == std::vector<std::string>{"x", "z"});
        CHECK(selected.col(0)[0] == 1.0);
        CHECK(selected.col(1)[0] == 7.0);
    }

    SECTION("invalid name throws") {
        CHECK_THROWS(named.select_columns({"x", "invalid"}));
    }
}

TEST_CASE("NamedDataset::save") {
    std::vector<double> x_data = {1.0, 2.0, 3.0};
    std::vector<double> y_data = {4.0, 5.0, 6.0};
    std::vector<double> z_data = {7.0, 8.0, 9.0};
    Dataset d({x_data, y_data, z_data});
    NamedDataset named(std::move(d), {"x", "y", "z"});

    SECTION("save to file") {
        io::File path("/tmp/named_dataset_test.dat");
        named.save(path, "# Test header");
        CHECK(path.exists());
        std::remove(path.str().c_str());
    }
}