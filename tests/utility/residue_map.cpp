#include <catch2/catch_test_macros.hpp>

#include <residue/ResidueMap.h>

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
        CHECK(residue_map.get_map().empty());
    }

    SECTION("unordered_map<AtomKey, int>&") {
        std::unordered_map<AtomKey, int> map;
        ResidueMap residue_map(map);
        CHECK(residue_map.get_map() == map);
    }
}

TEST_CASE("ResidueMap::get") {
    SECTION("AtomKey&") {
        ResidueMap residue_map;
        AtomKey key("name", constants::atom_t::Al);
        CHECK_THROWS(residue_map.get(key));

        residue_map.insert(key, 2);
        CHECK(residue_map.get(key) == 2);
    }

    SECTION("string&, string&") {
        ResidueMap residue_map;
        std::string name = "name";
        auto symbol = constants::atom_t::Al;
        CHECK_THROWS(residue_map.get(name, symbol));

        residue_map.insert(name, symbol, 2);
        CHECK(residue_map.get(name, symbol) == 2);
    }
}

TEST_CASE("ResidueMap::insert") {
    SECTION("AtomKey&, int") {
        ResidueMap residue_map;
        AtomKey key("name", constants::atom_t::Al);
        residue_map.insert(key, 2);
        CHECK(residue_map.get(key) == 2);
    }

    SECTION("string&, string&, int") {
        ResidueMap residue_map;
        std::string name = "name";
        auto symbol = constants::atom_t::Al;
        residue_map.insert(name, symbol, 2);
        CHECK(residue_map.get(name, symbol) == 2);
    }
}

TEST_CASE("ResidueMap::iterators") {
    SECTION("begin") {
        ResidueMap residue_map;
        CHECK(residue_map.begin() == residue_map.get_map().begin());
    }

    SECTION("end") {
        ResidueMap residue_map;
        CHECK(residue_map.end() == residue_map.get_map().end());
    }
}