#include <catch2/catch_test_macros.hpp>

#include <constants/Constants.h>
#include <io/File.h>
#include <io/detail/structure/CIFReader.h>
#include <residue/ResidueStorage.h>
#include <residue/detail/Residue.h>
#include <residue/detail/ResidueMap.h>
#include <utility/Curl.h>

#include <map>
#include <string>
#include <utility>

using namespace ausaxs;

TEST_CASE("ResidueParser: CCD definitions") {
    if (auto folder = io::Folder("temp/residues"); !folder.exists()) {
        folder.create();
    }

    // hydrogens per atom: {as the CCD defines them, in the pH 7 state the storage holds them}
    const std::map<std::string, std::map<std::string, std::pair<int, int>>> expected = {
        {"ASP", {{"N", {1, 1}}, {"CA", {1, 1}}, {"C", {0, 0}}, {"O", {0, 0}}, {"OXT", {1, 0}}, {"CB", {2, 2}}, {"CG", {0, 0}},
                 {"OD1", {0, 0}}, {"OD2", {1, 0}}}},
        {"HIS", {{"N", {1, 1}}, {"CA", {1, 1}}, {"C", {0, 0}}, {"O", {0, 0}}, {"OXT", {1, 0}}, {"CB", {2, 2}}, {"CG", {0, 0}},
                 {"ND1", {1, 0}}, {"CD2", {1, 1}}, {"CE1", {1, 1}}, {"NE2", {1, 1}}}}
    };
    auto element = [] (const std::string& atom) {return constants::symbols::parse_element_string(std::string(1, atom[0]));};

    // insert all before the first lookup, which adds the master definitions of the rest
    residue::ResidueStorage storage;
    for (const auto& [name, atoms] : expected) {
        io::File file("temp/residues/" + name + ".cif");
        if (!file.exists()) {
            REQUIRE(curl::download("files.rcsb.org/ligands/view/" + name + ".cif", file.str()));
        }

        auto parsed = io::detail::cif::read_residue(file).front();
        auto ccd = parsed.to_map();
        for (const auto& [atom, hydrogens] : atoms) {
            INFO("CCD " + name + " " + atom);
            CHECK(hydrogens.first == ccd.get(atom, element(atom)));
        }
        storage.insert_ccd(name, parsed.to_map());
    }

    for (const auto& [name, atoms] : expected) {
        auto& ph7 = storage.get(name);
        for (const auto& [atom, hydrogens] : atoms) {
            INFO("pH 7 " + name + " " + atom);
            CHECK(hydrogens.second == ph7.get(atom, element(atom)));
        }
    }
}
