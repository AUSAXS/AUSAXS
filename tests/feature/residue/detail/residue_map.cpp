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

TEST_CASE("ResidueMap::get: CHARMM C-terminal oxygens") {
    // CHARMM names the two C-terminal oxygens OT1/OT2, where the residue definitions carry the GROMACS names OC1/OC2
    ResidueMap residue_map({
        {AtomKey("OC1", constants::atom_t::O), 1},
        {AtomKey("OC2", constants::atom_t::O), 0},
        {AtomKey("O", constants::atom_t::O), 0}
    });
    CHECK(residue_map.get("OT1", constants::atom_t::O) == 1);
    CHECK(residue_map.get("OT2", constants::atom_t::O) == 0);
    CHECK(residue_map.get_atomic_group("OT1", constants::atom_t::O) == residue_map.get_atomic_group("OC1", constants::atom_t::O));
    CHECK(residue_map.get_atomic_group("OT2", constants::atom_t::O) == residue_map.get_atomic_group("OC2", constants::atom_t::O));

    // a name the definitions do carry is never replaced by its alias
    residue_map.get_backing_map().emplace(AtomKey("OT1", constants::atom_t::O), 0);
    CHECK(residue_map.get("OT1", constants::atom_t::O) == 0);
}
