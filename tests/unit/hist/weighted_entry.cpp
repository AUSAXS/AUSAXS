#include <catch2/catch_test_macros.hpp>

#include <hist/distribution/detail/WeightedEntry.h>

using namespace ausaxs::hist::detail;

TEST_CASE("WeightedEntry::WeightedEntry") {
    SECTION("default constructor") {
        WeightedEntry entry;
        CHECK(entry.value == 0);
        CHECK(entry.count == 0);
        CHECK(entry.bin_center == 0);
    }

    SECTION("parameterized constructor") {
        WeightedEntry entry(10, 5, 2.5);
        CHECK(entry.value == 10);
        CHECK(entry.count == 5);
        CHECK(entry.bin_center == 2.5);
    }
}

TEST_CASE("WeightedEntry::add") {
    SECTION("N=1") {
        WeightedEntry entry;
        entry.add<1>(5.0, 10);
        CHECK(entry.value == 10);
        CHECK(entry.count == 1);
        CHECK(entry.bin_center == 5.0);
    }

    SECTION("N=2") {
        WeightedEntry entry;
        entry.add<2>(5.0, 10);
        CHECK(entry.value == 20);
        CHECK(entry.count == 2);
        CHECK(entry.bin_center == 10.0);
    }

    SECTION("multiple additions") {
        WeightedEntry entry;
        entry.add<1>(2.0, 5);
        entry.add<1>(4.0, 10);
        CHECK(entry.value == 15);
        CHECK(entry.count == 2);
        CHECK(entry.bin_center == 6.0);
    }
}

TEST_CASE("WeightedEntry::operator+") {
    WeightedEntry entry1(10, 5, 2.5);
    WeightedEntry entry2(20, 3, 1.5);
    
    auto result = entry1 + entry2;
    CHECK(result.value == 30);
    CHECK(result.count == 8);
    CHECK(result.bin_center == 4.0);
}

TEST_CASE("WeightedEntry::operator+=") {
    WeightedEntry entry1(10, 5, 2.5);
    WeightedEntry entry2(20, 3, 1.5);
    
    entry1 += entry2;
    CHECK(entry1.value == 30);
    CHECK(entry1.count == 8);
    CHECK(entry1.bin_center == 4.0);
}

TEST_CASE("WeightedEntry::operator-") {
    WeightedEntry entry1(30, 8, 4.0);
    WeightedEntry entry2(20, 3, 1.5);
    
    auto result = entry1 - entry2;
    CHECK(result.value == 10);
    CHECK(result.count == 5);
    CHECK(result.bin_center == 2.5);
}

TEST_CASE("WeightedEntry::operator-=") {
    WeightedEntry entry1(30, 8, 4.0);
    WeightedEntry entry2(20, 3, 1.5);
    
    entry1 -= entry2;
    CHECK(entry1.value == 10);
    CHECK(entry1.count == 5);
    CHECK(entry1.bin_center == 2.5);
}

TEST_CASE("WeightedEntry::operator==") {
    WeightedEntry entry(10, 5, 2.5);
    CHECK(entry == 10);
    CHECK_FALSE(entry == 20);
}

TEST_CASE("WeightedEntry::operator*") {
    SECTION("entry * factor") {
        WeightedEntry entry(10, 5, 2.5);
        auto result = entry * 2.0;
        CHECK(result.value == 20);
        CHECK(result.count == 10);
        CHECK(result.bin_center == 5.0);
    }

    SECTION("factor * entry") {
        WeightedEntry entry(10, 5, 2.5);
        auto result = 2.0 * entry;
        CHECK(result.value == 20);
        CHECK(result.count == 10);
        CHECK(result.bin_center == 5.0);
    }
}
