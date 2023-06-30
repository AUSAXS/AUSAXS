#include <catch2/catch_test_macros.hpp>

#include <utility/Table.h>

using namespace table;

TEST_CASE("Table::Table") {
    SECTION("default") {
        Table table;
        CHECK(table.data.empty());
        CHECK(table.N == 0);
        CHECK(table.M == 0);
    }

    SECTION("N, M") {
        Table table(3, 4);
        CHECK_FALSE(table.data.empty());
        CHECK(table.N == 3);
        CHECK(table.M == 4);
    }
}

TEST_CASE("Table::index") {
    Table table(3, 4);
    CHECK(table.index(0, 0) == 0);
    CHECK(table.index(2, 3) == 0);
}