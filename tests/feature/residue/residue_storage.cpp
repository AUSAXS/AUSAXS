#include <catch2/catch_test_macros.hpp>

#include <residue/ResidueStorage.h>
#include <residue/detail/ResidueMap.h>
#include <settings/MoleculeSettings.h>

#include <algorithm>
#include <initializer_list>
#include <string>

using namespace ausaxs;
using namespace residue::detail;

TEST_CASE("ResidueStorage: atom name aliases") {
    auto C = constants::atom_t::C;
    auto O = constants::atom_t::O;

    residue::ResidueStorage storage;
    storage.insert("ILE", ResidueMap({{AtomKey("O", O), 0}, {AtomKey("OXT", O), 1}, {AtomKey("CD1", C), 3}}));
    storage.insert("TST", ResidueMap({{AtomKey("O", O), 0}, {AtomKey("OXT", O), 1}, {AtomKey("CD1", C), 3}, {AtomKey("CA", C), 1}, {AtomKey("OT1", O), 0}}));

    SECTION("C-terminal oxygens") {
        for (const auto* name : {"OC1", "OT1"}) {
            CHECK(storage.get_hydrogens("ILE", name, O) == 1);
            CHECK(storage.get_atomic_group("ILE", name, O) == storage.get_atomic_group("ILE", "OXT", O));
        }
        for (const auto* name : {"OC2", "OT2"}) {
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


TEST_CASE("ResidueStorage: force-field protonation-state names") {
    settings::molecule::allow_unknown_atoms = false;
    auto N = constants::atom_t::N;
    auto O = constants::atom_t::O;
    auto S = constants::atom_t::S;
    residue::ResidueStorage storage;
    storage.insert("HSE", ResidueMap({{AtomKey("OG", O), 1}})); // stands in for the CCD definition, homoserine

    SECTION("the CCD definition by default") {
        settings::molecule::force_field_residue_names = false;
        CHECK(storage.get_hydrogens("HSE", "OG", O) == 1);
        CHECK_THROWS(storage.get_atomic_group("HSE", "NE2", N));
    }

    SECTION("the force-field state if enabled") {
        settings::molecule::force_field_residue_names = true;
        CHECK(storage.get_hydrogens("HSE", "NE2", N) == 1);
        CHECK_THROWS(storage.get_atomic_group("HSE", "OG", O));

        for (const auto* res : {"HID", "HSD"}) {CHECK(storage.get_hydrogens(res, "ND1", N) == 1); CHECK(storage.get_hydrogens(res, "NE2", N) == 0);}
        for (const auto* res : {"HIE", "HSE"}) {CHECK(storage.get_hydrogens(res, "ND1", N) == 0); CHECK(storage.get_hydrogens(res, "NE2", N) == 1);}
        for (const auto* res : {"HIP", "HSP"}) {CHECK(storage.get_hydrogens(res, "ND1", N) == 1); CHECK(storage.get_hydrogens(res, "NE2", N) == 1);}
        for (const auto* res : {"CYX", "CYM"}) {CHECK(storage.get_hydrogens(res, "SG", S) == 0);}
        CHECK(storage.get_hydrogens("ASH", "OD1", O) == 0); CHECK(storage.get_hydrogens("ASH", "OD2", O) == 1);
        CHECK(storage.get_hydrogens("GLH", "OE1", O) == 0); CHECK(storage.get_hydrogens("GLH", "OE2", O) == 1);
        CHECK(storage.get_hydrogens("LYN", "NZ", N) == 2);

        // every other atom is that of the standard residue
        auto check = [&] (const std::string& variant, const std::string& base, std::initializer_list<std::string> titratable) {
            for (const auto& [key, hydrogens] : storage.get(base).get_backing_map()) {
                if (std::ranges::find(titratable, key.name) != titratable.end()) {continue;}
                CHECK(storage.get_hydrogens(variant, key.name, key.atom) == hydrogens);
            }
        };
        for (const auto* res : {"HID", "HSD", "HIE", "HSE", "HIP", "HSP"}) {check(res, "HIS", {"nd1", "ne2"});}
        for (const auto* res : {"CYX", "CYM"}) {check(res, "CYS", {"sg"});}
        check("ASH", "ASP", {"od1", "od2"});
        check("GLH", "GLU", {"oe1", "oe2"});
        check("LYN", "LYS", {"nz"});
    }
}
