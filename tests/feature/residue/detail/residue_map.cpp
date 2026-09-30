#include <catch2/catch_test_macros.hpp>

#include <residue/detail/ResidueMap.h>

using namespace ausaxs;
using namespace residue::detail;

TEST_CASE("AtomKey::AtomKey") {
    SECTION("string&, string&") {
        std::string name = "name";
        auto symbol = constants::atom_t::Al;
        AtomKey key(name, symbol);
        CHECK(key.name == name);
        CHECK(key.atom == symbol);
    }
}

TEST_CASE("AtomKey::operator==") {
    SECTION("const AtomKey&") {
        AtomKey key1("name", constants::atom_t::Al);
        AtomKey key2("name", constants::atom_t::Al);
        CHECK(key1 == key2);

        AtomKey key3("name", constants::atom_t::H);
        CHECK(key1 == key3);

        AtomKey key4("name2", constants::atom_t::Al);
        CHECK_FALSE(key1 == key4);
    }
}

TEST_CASE("ResidueMap::ResidueMap") {
    SECTION("default") {
        ResidueMap residue_map;
        CHECK(residue_map.get_backing_map().empty());
    }

    SECTION("unordered_map<AtomKey, int>&") {
        std::unordered_map<AtomKey, int> map;
        ResidueMap residue_map(map);
        CHECK(residue_map.get_backing_map() == map);
    }
}