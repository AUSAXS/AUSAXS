#include <catch2/catch_test_macros.hpp>

#include <table/Table.h>

using namespace table;

TEST_CASE("Table::Table") {
    SECTION("default") {
        Table table;
        CHECK(table.size_x() == 0);
        CHECK(table.size_y() == 0);
    }

    SECTION("N, M") {
        Table table(3, 4);
        CHECK(table.size_x() == 3);
        CHECK(table.size_y() == 4);
    }
}

TEST_CASE("Table::index") {
    Table table(3, 4);
    CHECK(table.index(0, 0) == 0);
    CHECK(table.index(2, 3) == 0);
}