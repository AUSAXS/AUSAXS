#include <catch2/catch_test_macros.hpp>

#include <constants/Constants.h>
#include <residue/ResidueStorage.h>
#include <residue/detail/ResidueMap.h>

using namespace ausaxs;
using namespace residue::detail;

TEST_CASE("ResidueStorage: atom name aliases") {
    auto C = constants::atom_t::C;
    auto O = constants::atom_t::O;

    residue::ResidueStorage storage;
    storage.insert("ILE", ResidueMap({{AtomKey("O", O), 0}, {AtomKey("OXT", O), 1}, {AtomKey("CD1", C), 3}}));
    storage.insert("TST", ResidueMap({{AtomKey("O", O), 0}, {AtomKey("OXT", O), 1}, {AtomKey("CD1", C), 3}, {AtomKey("CA", C), 1}, {AtomKey("OT1", O), 0}}));

    SECTION("C-terminal oxygens") {
        for (auto name : {"OC1", "OT1"}) {
            CHECK(storage.get_hydrogens("ILE", name, O) == 1);
            CHECK(storage.get_atomic_group("ILE", name, O) == storage.get_atomic_group("ILE", "OXT", O));
        }
        for (auto name : {"OC2", "OT2"}) {
            CHECK(storage.get_hydrogens("ILE", name, O) == 0);
            CHECK(storage.get_atomic_group("ILE", name, O) == storage.get_atomic_group("ILE", "O", O));
        }
    }

    SECTION("residue-specific") {
        CHECK(storage.get_hydrogens("ILE", "CD", C) == 3);
        CHECK(storage.get_hydrogens("TST", "CD", C) == 2); // unknown name: the carbon average
    }

    SECTION("a name the definition carries is not an alias") {
        CHECK(storage.get_hydrogens("TST", "OT1", O) == 0);
    }
}
