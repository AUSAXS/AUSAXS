#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Water.h>
#include <data/record/Record.h>
#include <utility/Console.h>
#include <settings/All.h>
#include <io/CIFReader.h>
#include <residue/ResidueParser.h>
#include <constants/Constants.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace ausaxs;
using namespace data;
using namespace data::record;

TEST_CASE("CIFReader::read") {
    settings::molecule::center = false;
    settings::general::verbose = false;

    io::File path("temp/io/temp.cif");
    path.create();

    std::ofstream cif_file(path);
    cif_file << "loop_" << std::endl;
    cif_file << "_atom_site.group_PDB" << std::endl;
    cif_file << "_atom_site.id" << std::endl;
    cif_file << "_atom_site.type_symbol" << std::endl;
    cif_file << "_atom_site.label_atom_id" << std::endl;
    cif_file << "_atom_site.label_alt_id" << std::endl;
    cif_file << "_atom_site.label_comp_id" << std::endl;
    cif_file << "_atom_site.label_asym_id" << std::endl;
    cif_file << "_atom_site.label_entity_id" << std::endl;
    cif_file << "_atom_site.label_seq_id" << std::endl;
    cif_file << "_atom_site.pdbx_PDB_ins_code" << std::endl;
    cif_file << "_atom_site.Cartn_x" << std::endl;
    cif_file << "_atom_site.Cartn_y" << std::endl;
    cif_file << "_atom_site.Cartn_z" << std::endl;
    cif_file << "_atom_site.occupancy" << std::endl;
    cif_file << "_atom_site.B_iso_or_equiv" << std::endl;
    cif_file << "_atom_site.pdbx_formal_charge" << std::endl;
    cif_file << "_atom_site.auth_seq_id" << std::endl;
    cif_file << "_atom_site.auth_comp_id" << std::endl;
    cif_file << "_atom_site.auth_asym_id" << std::endl;
    cif_file << "_atom_site.auth_atom_id" << std::endl;
    cif_file << "_atom_site.pdbx_PDB_model_num" << std::endl;
    cif_file << "_atom_site.calc_flag" << std::endl;
    cif_file << "ATOM   1    N N   . SER A 1 1   ? -32.928 -5.043  -33.904 0.100 45.216  0 1   SER AAA N   1 ?" << std::endl;
    cif_file << "ATOM   2    C CA  . SER A 1 1   ? -32.489 -6.122  -32.994 0.200 46.910  0 1   SER AAA CA  1 ?" << std::endl;
    cif_file << "ATOM   3    C C   . SER A 1 1   ? -30.974 -6.329  -33.145 0.300 46.121  0 1   SER AAA C   1 ?" << std::endl;
    cif_file << "ATOM   4    O O   . SER A 1 1   ? -30.318 -5.462  -33.746 0.400 50.538  0 1   SER AAA O   1 ?" << std::endl;
    cif_file << "#" << std::endl;
    cif_file.close();

    // check CIF io
    data::Molecule protein(path);
    auto atoms = protein.get_atoms();
    REQUIRE(atoms.size() == 4);
    REQUIRE(atoms[0].get_serial() == 1);
    REQUIRE(atoms[0].get_coordinates().x() == -32.928);
    REQUIRE(atoms[0].get_coordinates().y() == -5.043);
    REQUIRE(atoms[0].get_coordinates().z() == -33.904);
    REQUIRE(atoms[0].get_occupancy() == 0.1);
    REQUIRE(atoms[0].get_element() == constants::atom_t::N);
    REQUIRE(atoms[0].get_residue_name() == "SER");
    REQUIRE(atoms[0].get_residue_sequence_number() == 1);
    REQUIRE(atoms[0].get_chainID() == 'A');
    REQUIRE(atoms[0].get_temperature_factor() == 45.216);

    REQUIRE(atoms[1].get_serial() == 2);
    REQUIRE(atoms[1].get_coordinates().x() == -32.489);
    REQUIRE(atoms[1].get_coordinates().y() == -6.122);
    REQUIRE(atoms[1].get_coordinates().z() == -32.994);
    REQUIRE(atoms[1].get_occupancy() == 0.2);
    REQUIRE(atoms[1].get_element() == constants::atom_t::C);
    REQUIRE(atoms[1].get_residue_name() == "SER");
    REQUIRE(atoms[1].get_residue_sequence_number() == 1);
    REQUIRE(atoms[1].get_chainID() == 'A');
    REQUIRE(atoms[1].get_temperature_factor() == 46.910);

    REQUIRE(atoms[2].get_serial() == 3);
    REQUIRE(atoms[2].get_coordinates().x() == -30.974);
    REQUIRE(atoms[2].get_coordinates().y() == -6.329);
    REQUIRE(atoms[2].get_coordinates().z() == -33.145);
    REQUIRE(atoms[2].get_occupancy() == 0.3);
    REQUIRE(atoms[2].get_element() == constants::atom_t::C);
    REQUIRE(atoms[2].get_residue_name() == "SER");
    REQUIRE(atoms[2].get_residue_sequence_number() == 1);
    REQUIRE(atoms[2].get_chainID() == 'A');
    REQUIRE(atoms[2].get_temperature_factor() == 46.121);

    REQUIRE(atoms[3].get_serial() == 4);
    REQUIRE(atoms[3].get_coordinates().x() == -30.318);
    REQUIRE(atoms[3].get_coordinates().y() == -5.462);
    REQUIRE(atoms[3].get_coordinates().z() == -33.746);
    REQUIRE(atoms[3].get_occupancy() == 0.4);
    REQUIRE(atoms[3].get_element() == constants::atom_t::O);
    REQUIRE(atoms[3].get_residue_name() == "SER");
    REQUIRE(atoms[3].get_residue_sequence_number() == 1);
    REQUIRE(atoms[3].get_chainID() == 'A');
    REQUIRE(atoms[3].get_temperature_factor() == 50.538);
}

TEST_CASE("CIFReader: uses file residues") {
    settings::general::verbose = false;

    data::Molecule cif("tests/files/3sba.cif");
    auto residues = io::detail::CIFReader::read_residue("tests/files/3sba.cif");

    for (const auto& residue : residues) {
        auto& loaded = constants::hydrogen_atoms::residues.get(residue.get_name());
        auto expected = residue.to_map();
        for (const auto& [atom, num] : expected) {
            REQUIRE(loaded.get(atom) == num);
        }
    }
}

TEST_CASE("CIFReader: file residues agrees with PDB") {
    settings::general::verbose = false;

    // loading the PDB file first to load default ResidueStorage data
    data::Molecule pdb("tests/files/3sba.pdb");
    pdb.add_implicit_hydrogens(); // ensure hydrogen data is initialized

    auto residues = io::detail::CIFReader::read_residue("tests/files/3sba.cif");
    for (const auto& residue : residues) {
        auto& loaded = constants::hydrogen_atoms::residues.get(residue.get_name());
        auto expected = residue.to_map();
        for (const auto& [atom, num] : expected) {
            if (loaded.get(atom) != num) {
                std::cout << "Residue: " << residue.get_name() << " Atom: " << constants::symbols::to_string(atom.atom) << " Expected: " << num << " Loaded: " << loaded.get(atom) << std::endl;
            }
            REQUIRE(loaded.get(atom) == num);
        }
    }
}

TEST_CASE("CIFReader: compare with PDB", "[files]") {
    settings::general::verbose = false;

    data::Molecule cif1("tests/files/3sba.cif");
    data::Molecule cif2("tests/files/7xb3.cif");
    data::Molecule cif3("tests/files/168l.cif");

    data::Molecule pdb1("tests/files/3sba.pdb");
    data::Molecule pdb2("tests/files/7xb3.pdb");
    data::Molecule pdb3("tests/files/168l.pdb");

    auto compare_atoms = [] (std::vector<Atom>&& a1, std::vector<Atom>&& a2) {
        REQUIRE(a1.size() == a2.size());
        int chain = 0;
        char chainID = a2[0].get_chainID();
        for (unsigned int i = 0; i < a1.size(); i++) {
            // if the chainID changes, the PDB file may or may not have a TER record, thus shifting all following serials by one
            if (a2[i].get_chainID() != chainID) {
                // check if the PDB atom is shifted relative to the CIF file
                if (a2[i].get_serial()-chain != a1[i].get_serial()) {
                    // if so, increment the chain shift (this is typically the case for new ATOM chains)
                    chain++;
                }
                // otherwise, just update the chainID (this is typically the case for HETATM)
                chainID = a2[i].get_chainID();
            }
            a2[i].serial -= chain; // account for the chain shift

            REQUIRE(a1[i].get_serial() == a2[i].get_serial());
            REQUIRE_THAT(a1[i].get_coordinates().x(), Catch::Matchers::WithinAbsMatcher(a2[i].get_coordinates().x(), 1e-2));
            REQUIRE_THAT(a1[i].get_coordinates().y(), Catch::Matchers::WithinAbsMatcher(a2[i].get_coordinates().y(), 1e-2));
            REQUIRE_THAT(a1[i].get_coordinates().z(), Catch::Matchers::WithinAbsMatcher(a2[i].get_coordinates().z(), 1e-2));
            REQUIRE_THAT(a1[i].get_occupancy(), Catch::Matchers::WithinAbs(a2[i].get_occupancy(), 1e-2));
            REQUIRE_THAT(a1[i].get_temperature_factor(), Catch::Matchers::WithinAbs(a2[i].get_temperature_factor(), 1e-2));
            REQUIRE(a1[i].get_element() == a2[i].get_element());
            REQUIRE(a1[i].get_residue_name() == a2[i].get_residue_name());
            REQUIRE(a1[i].get_residue_sequence_number() == a2[i].get_residue_sequence_number());
            REQUIRE(a1[i].get_chainID() == a2[i].get_chainID());
            if (!a1[i].get_alternate_location().starts_with('.')) {REQUIRE(a1[i].get_alternate_location() == a2[i].get_alternate_location());}
        }
    };

    compare_atoms(cif1.get_atoms(), pdb1.get_atoms());
    compare_atoms(cif2.get_atoms(), pdb2.get_atoms());
    compare_atoms(cif3.get_atoms(), pdb3.get_atoms());
}