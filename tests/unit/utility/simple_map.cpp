#include <catch2/catch_test_macros.hpp>

#include <utility/SimpleMap.h>

using namespace ausaxs;
using namespace saxs::detail;

TEST_CASE("SimpleMap::SimpleMap") {
    SECTION("default") {
        SimpleMap<int> map;
        CHECK(map.get_map().empty());
    }

    SECTION("unordered_map<string, V>&") {
        std::unordered_map<std::string, int> map = {{"key1", 1}, {"key2", 2}};
        SimpleMap<int> simple_map(map);
        CHECK(simple_map.get_map() == map);
    }
}

TEST_CASE("SimpleMap::get") {
    std::unordered_map<std::string, int> map = {{"key1", 1}, {"key2", 2}, {"KeY3", 3}};
    SimpleMap<int> simple_map(map);
    CHECK(simple_map.get("key1") == 1);
    CHECK(simple_map.get("key2") == 2);
    CHECK(simple_map.get("key3") == 3);
    CHECK_THROWS(simple_map.get("key4"));
}

TEST_CASE("SimpleMap::insert") {
    SimpleMap<int> simple_map;
    simple_map.insert("key1", 1);
    simple_map.insert("key2", 2);
    simple_map.insert("KeY3", 3);
    CHECK(simple_map.get("key1") == 1);
    CHECK(simple_map.get("key2") == 2);
    CHECK(simple_map.get("key3") == 3);
}

TEST_CASE("SimpleMap::contains") {
    std::unordered_map<std::string, int> map = {{"key1", 1}, {"key2", 2}, {"KeY3", 3}};
    SimpleMap<int> simple_map(map);
    CHECK(simple_map.contains("key1"));
    CHECK(simple_map.contains("key2"));
    CHECK(simple_map.contains("key3"));
    CHECK_FALSE(simple_map.contains("key4"));
}

TEST_CASE("SimpleMap::iterators") {
    std::unordered_map<std::string, int> map = {{"key1", 1}, {"key2", 2}, {"KeY3", 3}};
    SimpleMap<int> simple_map(map);

    int i = 0;
    for (auto it = simple_map.begin(); it != simple_map.end(); ++it) {++i;}
    CHECK(i == 3);
}