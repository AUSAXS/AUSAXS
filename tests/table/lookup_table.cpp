#include <catch2/catch_test_macros.hpp>

#include <table/LookupTable.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;

TEST_CASE("LookupTable::LookupTable") {
    SECTION("default") {
        table::LookupTable<int, int> table;
        CHECK(table.is_empty());
    }

    SECTION("std::vector<T>&, std::vector<Q>&") {
        std::vector<int> rows = {1, 2, 3};
        std::vector<int> cols = {4, 5, 6};
        table::LookupTable<int, int> table(rows, cols);
        CHECK(table.is_empty() == false);
    }
}

TEST_CASE("LookupTable::initialize") {
    SECTION("basic") {
        std::vector<int> rows = {1, 2, 3};
        std::vector<int> cols = {4, 5, 6};
        table::LookupTable<int, int> table;
        table.initialize(rows, cols);
        CHECK(table.is_empty() == false);
    }

    SECTION("strings") {
        std::vector<std::string> rows = {"a", "b", "c"};
        std::vector<std::string> cols = {"d", "e", "f"};
        table::LookupTable<std::string, std::string> table;
        table.initialize(rows, cols);
        CHECK(table.is_empty() == false);

        CHECK(table.lookup("a", "d") == 0);
        CHECK(table.lookup("a", "e") == 0);
    }
}

TEST_CASE("LookupTable::is_empty") {
    SECTION("empty") {
        table::LookupTable<int, int> table;
        CHECK(table.is_empty());
    }

    SECTION("not empty") {
        std::vector<int> rows = {1, 2, 3};
        std::vector<int> cols = {4, 5, 6};
        table::LookupTable<int, int> table(rows, cols);
        CHECK(table.is_empty() == false);
    }
}

TEST_CASE("LookupTable::assign_index") {
    std::vector<std::string> rows = {"a", "b", "c"};
    std::vector<std::string> cols = {"d", "e", "f"};

    table::LookupTable<std::string, std::string> table(rows, cols);
    CHECK(table.lookup_index(0, 1) == 0);
    CHECK(table.lookup_index(1, 2) == 0);

    table.assign_index(0, 1, 1);
    CHECK(table.lookup_index(0, 1) == 1);

    table.assign_index(1, 2, 2);
    CHECK(table.lookup_index(1, 2) == 2);
}

TEST_CASE("LookupTable::assign") {
    std::vector<std::string> rows = {"a", "b", "c"};
    std::vector<std::string> cols = {"d", "e", "f"};

    table::LookupTable<std::string, std::string> table(rows, cols);
    CHECK(table.lookup("a", "d") == 0);
    CHECK(table.lookup("a", "e") == 0);

    table.assign("a", "d", 1);
    CHECK(table.lookup("a", "d") == 1);
    CHECK(table.lookup_index(0, 0) == 1);

    table.assign("c", "f", 2);
    CHECK(table.lookup("c", "f") == 2);
    CHECK(table.lookup_index(2, 2) == 2);
}


TEST_CASE("lookup_tables") {
    settings::axes::qmax = 1.001;
    SECTION("integer table") {
        std::vector<unsigned int> rows = {1, 2, 3};
        std::vector<unsigned int> cols = {1, 2, 3};
        table::LookupTable<unsigned int, unsigned int> table(rows, cols);

        for (unsigned int r : rows) {
            for (unsigned int c : cols) {
                table.assign(r, c, r*c);
            }
        }

        CHECK(table.lookup(1, 1) == 1);
        CHECK(table.lookup(1, 2) == 2);
        CHECK(table.lookup(1, 3) == 3);
        CHECK(table.lookup(2, 1) == 2);
        CHECK(table.lookup(2, 2) == 4);
        CHECK(table.lookup(2, 3) == 6);
        CHECK(table.lookup(3, 1) == 3);
        CHECK(table.lookup(3, 2) == 6);
        CHECK(table.lookup(3, 3) == 9);
    }

    SECTION("double table") {
        std::vector<double> rows = {0.35, -2.5, 10.2};
        std::vector<double> cols = {-3.1, -5.05, 6.7};
        table::LookupTable<double, double> table(rows, cols);

        for (double r : rows) {
            for (double c : cols) {
                table.assign(r, c, r+c);
            }
        }

        CHECK(table.lookup(rows[0], cols[0]) == rows[0] + cols[0]);
        CHECK(table.lookup(rows[0], cols[1]) == rows[0] + cols[1]);
        CHECK(table.lookup(rows[0], cols[2]) == rows[0] + cols[2]);
        CHECK(table.lookup(rows[1], cols[0]) == rows[1] + cols[0]);
        CHECK(table.lookup(rows[1], cols[1]) == rows[1] + cols[1]);
        CHECK(table.lookup(rows[1], cols[2]) == rows[1] + cols[2]);
        CHECK(table.lookup(rows[2], cols[0]) == rows[2] + cols[0]);
        CHECK(table.lookup(rows[2], cols[1]) == rows[2] + cols[1]);
        CHECK(table.lookup(rows[2], cols[2]) == rows[2] + cols[2]);        
    }
}