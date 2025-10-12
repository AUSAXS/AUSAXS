#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <constants/Constants.h>
#include <utility/Curl.h>

#include <map>
#include <string>

using namespace ausaxs;

/**
 * @brief The old version of the Constants.h file. This is used for testing the new version.
 */
namespace hydrogen_atoms {
    namespace glycine {
        constexpr int N = 1;
        constexpr int CA = 2;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}};
    }
    namespace alanine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 3;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}};
    }
    namespace valine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 1;
        constexpr int CG1 = 3;
        constexpr int CG2 = 3;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG1", CG1}, {"CG2", CG2}};
    }
    namespace leucine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 1;
        constexpr int CD1 = 3;
        constexpr int CD2 = 3;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, {"CD2", CD2}};
    }
    namespace isoleucine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 1;
        constexpr int CG2 = 3;
        constexpr int CG1 = 2;
        constexpr int CD1 = 3;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG2", CG2}, {"CG1", CG1}, {"CD1", CD1}};
    }
    namespace phenylalanine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 0;
        constexpr int CD1 = 1;
        constexpr int CD2 = 1;
        constexpr int CE1 = 1;
        constexpr int CE2 = 1;
        constexpr int CZ = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, 
            {"CD2", CD2}, {"CE1", CD1}, {"CE2", CD2}, {"CZ", CZ}};
    }
    namespace tyrosine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 0;
        constexpr int CD1 = 1;
        constexpr int CD2 = 1;
        constexpr int CE1 = 1;
        constexpr int CE2 = 1;
        constexpr int CZ = 0;
        constexpr int OH = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, 
            {"CD2", CD2}, {"CE1", CE1}, {"CE2", CE2}, {"CZ", CZ}, {"OH", OH}};
    }
    namespace tryptophan {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 0;
        constexpr int CD1 = 1;
        constexpr int CD2 = 0;
        constexpr int NE1 = 1;
        constexpr int CE2 = 0;
        constexpr int CE3 = 1;
        constexpr int CZ2 = 1;
        constexpr int CZ3 = 1;
        constexpr int CH2 = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD1", CD1}, 
            {"CD2", CD2}, {"NE1", NE1}, {"CE2", CE2}, {"CE3", CE3}, {"CZ2", CZ2}, {"CZ3", CZ3}, {"CH2", CH2}};
    }
    namespace aspartic_acid {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 0;
        constexpr int OD1 = 0;
        constexpr int OD2 = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"OD1", OD1}, 
            {"OD2", OD2}};
    }
    namespace glutamic_acid {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 2;
        constexpr int CD = 0;
        constexpr int OE1 = 0;
        constexpr int OE2 = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
            {"OE1", OE1}, {"OE2", OE2}};
    }
    namespace serine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int OG = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"OG", OG}};
    }
    namespace threonine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 1;
        constexpr int OG1 = 1;
        constexpr int CG2 = 3;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"OG1", OG1}, {"CG2", CG2}};
    }
    namespace asparagine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 0;
        constexpr int OD1 = 0;
        constexpr int ND2 = 2;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"OD1", OD1}, 
            {"ND2", ND2}};
    }
    namespace glutamine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 2;
        constexpr int CD = 0;
        constexpr int OE1 = 0;
        constexpr int NE2 = 2;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
            {"OE1", OE1}, {"NE2", NE2}};
    }
    namespace lysine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 2;
        constexpr int CD = 2;
        constexpr int CE = 2;
        constexpr int NZ = 3;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
            {"CE", CE}, {"NZ", NZ}};
    }
    namespace arginine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 2;
        constexpr int CD = 2;
        constexpr int NE = 1;
        constexpr int CZ = 0;
        constexpr int NH1 = 2;
        constexpr int NH2 = 2;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}, 
            {"NE", NE}, {"CZ", CZ}, {"NH1", NH1}, {"NH2", NH2}};
    }
    namespace histidine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 0;
        constexpr int ND1 = 1;
        constexpr int CD2 = 1;
        constexpr int CE1 = 1;
        constexpr int NE2 = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"ND1", ND1}, 
            {"CD2", CD2}, {"CE1", CE1}, {"NE2", NE2}};
    }
    namespace methionine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 2;
        constexpr int SD = 0;
        constexpr int CE = 3;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"SD", SD}, 
            {"CE", CE}};
    }
    namespace cysteine {
        constexpr int N = 1;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int SG = 1;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"SG", SG}};
    }
    namespace proline {
        constexpr int N = 0;
        constexpr int CA = 1;
        constexpr int C = 0;
        constexpr int O = 0;
        constexpr int OXT = 1;
        constexpr int CB = 2;
        constexpr int CG = 2;
        constexpr int CD = 2;
        const std::map<std::string, unsigned int> get = {{"N", N}, {"CA", CA}, {"C", C}, {"O", O}, {"OXT", OXT}, {"CB", CB}, {"CG", CG}, {"CD", CD}};
    }
    namespace myristic_acid {
        constexpr int C1 = 0;
        constexpr int O1 = 0;
        constexpr int O2 = 1;
        constexpr int C2 = 2;
        constexpr int C3 = 2;
        constexpr int C4 = 2;
        constexpr int C5 = 2;
        constexpr int C6 = 2;
        constexpr int C7 = 2;
        constexpr int C8 = 2;
        constexpr int C9 = 2;
        constexpr int C10 = 2;
        constexpr int C11 = 2;
        constexpr int C12 = 2;
        constexpr int C13 = 2;
        constexpr int C14 = 3;
        const std::map<std::string, unsigned int> get = {{"C1", C1}, {"O1", O1}, {"O2", O2}, {"C2", C2}, {"C3", C3}, {"C4", C4}, {"C5", C5}, {"C6", C6}, {"C7", C7},
            {"C8", C8}, {"C9", C9}, {"C10", C10}, {"C11", C11}, {"C12", C12}, {"C13", C13}, {"C14", C14}};
    }

    // get the number of hydrogen atoms attached to an atom of a specific acid. Example: get.at("GLY").at("CA") = 2
    const std::map<std::string, std::map<std::string, unsigned int>> get = {{"GLY", glycine::get}, {"ALA", alanine::get}, {"VAL", valine::get}, 
        {"LEU", leucine::get}, {"ILE", isoleucine::get}, {"PHE", phenylalanine::get}, {"TYR", tyrosine::get}, {"TRP", tryptophan::get}, 
        {"ASP", aspartic_acid::get}, {"GLU", glutamic_acid::get}, {"SER", serine::get}, {"THR", threonine::get}, {"ASN", asparagine::get}, 
        {"GLN", glutamine::get}, {"LYS", lysine::get}, {"ARG", arginine::get}, {"HIS", histidine::get}, {"MET", methionine::get}, 
        {"CYS", cysteine::get}, {"PRO", proline::get}, {"MYR", myristic_acid::get}};
}

TEST_CASE("ResidueParser: parse_single") {
    auto GLY = constants::hydrogen_atoms::residues.get("GLY");
    auto GLY2 = hydrogen_atoms::glycine::get;

    for (const auto& [atom, num] : GLY2) {
        REQUIRE(num == GLY.get(atom, constants::symbols::parse_element_string(std::string(1, atom[0]))));
    }
}

TEST_CASE("ResidueParser: parse_all") {
    for (const auto& [acid, atom_map] : hydrogen_atoms::get) {
        for (const auto& [atom, num_hydrogens] : atom_map) {
            SECTION(acid + " " + atom) {
                CHECK(hydrogen_atoms::get.at(acid).at(atom) == constants::hydrogen_atoms::residues.get(acid).get(atom, constants::symbols::parse_element_string(std::string(1, atom[0]))));
            }
        }
    }
}

#include <io/detail/CIFReader.h>
#include <settings/GeneralSettings.h>
TEST_CASE("ResidueParser: cif_reader_single") {
    if (auto folder = io::Folder("temp/residues"); !folder.exists()) {
        folder.create();
    }

    io::File gly("temp/residues/GLY.cif");
    if (!gly.exists()) {
        REQUIRE(ausaxs::curl::download("files.rcsb.org/ligands/view/GLY.cif", "temp/residues/GLY.cif"));
    }
    auto residue = io::detail::cif::read_residue(gly).front().to_map();
    for (const auto& [atom, num] : hydrogen_atoms::glycine::get) {
        REQUIRE(num == residue.get(atom, constants::symbols::parse_element_string(std::string(1, atom[0]))));
    }
}

TEST_CASE("ResidueParser: cif_reader_all") {
    if (auto folder = io::Folder("temp/residues"); !folder.exists()) {
        folder.create();
    }

    for (const auto& [acid, atom_map] : hydrogen_atoms::get) {
        io::File res("temp/residues/" + acid + ".cif");
        if (!res.exists()) {
            REQUIRE(ausaxs::curl::download("files.rcsb.org/ligands/view/" + acid + ".cif", res.str()));
        }

        auto residue = io::detail::cif::read_residue(res).front().to_map();
        for (const auto& [atom, num_hydrogens] : atom_map) {
            SECTION(acid + " " + atom) {
                CHECK(hydrogen_atoms::get.at(acid).at(atom) == residue.get(atom, constants::symbols::parse_element_string(std::string(1, atom[0]))));
            }
        }
    }
}