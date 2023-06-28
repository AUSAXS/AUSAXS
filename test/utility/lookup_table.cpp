#include <catch2/catch_test_macros.hpp>

#include <utility/LookupTable.h>

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