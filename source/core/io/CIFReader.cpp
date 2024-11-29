/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/CIFReader.h>
#include <io/File.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <utility/Console.h>
#include <residue/ResidueParser.h>
#include <constants/Constants.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/detail/AtomCollection.h>
#include <settings/All.h>

#include <fstream>
#include <unordered_map>

using namespace ausaxs;

io::detail::CIFReader::CIFReader(data::detail::AtomCollection* const file) : file(file) {}

io::detail::CIFReader::~CIFReader() = default;

using CIFRow = std::vector<std::string>;
struct CIFSection {
    std::vector<std::string> labels;
    std::vector<CIFRow> data;

    bool empty() const {return labels.empty() && data.empty();}

    std::unordered_map<std::string, int> get_label_map() const {
        std::unordered_map<std::string, int> map;
        for (size_t i = 0; i < labels.size(); i++) {
            map[labels[i]] = i;
        }
        return map;
    }
};

residue::detail::Residue parse_ion(const CIFSection& atom) {
    if (atom.data.size() != 1) {throw except::io_error("CIFReader::parse_ion: Invalid number of data entries in atom section");}
    auto labels = atom.get_label_map();

    if (!labels.contains("comp_id")) {
        throw except::io_error("CIFReader::parse_ion: Missing required label \"comp_id\" in \"_chem_comp\" section");
    } if (!labels.contains("charge")) {
        throw except::io_error("CIFReader::parse_ion: Missing required label \"charge\" in \"_chem_comp\" section");
    }

    std::string comp_id = atom.data[0][labels.at("comp_id")];
    if (!constants::symbols::detail::string_to_atomt_map.contains(comp_id)) {
        throw except::io_error("CIFReader::parse_ion: Unrecognized ion: \"" + comp_id + "\"");
    }

    int charge = std::stoi(atom.data[0][labels.at("charge")]);
    residue::detail::Residue residue(comp_id);
    residue.add_atom(comp_id, charge, constants::symbols::parse_element_string(comp_id));
    return residue;
}

std::vector<residue::detail::Residue> parse_residue(CIFSection& atom, CIFSection& bond) {
    if (atom.data.empty() || bond.data.empty()) {throw except::io_error("CIFReader::parse_residue: Empty data section");}
    auto atom_labels = atom.get_label_map();
    auto bond_labels = bond.get_label_map();

    // mandatory labels
    if (!atom_labels.contains("comp_id") || !atom_labels.contains("atom_id")) {
        throw except::io_error("CIFReader::parse_residue: Missing required labels in \"_chem_comp_atom\" section");
    }

    // optional label
    bool has_alt_atom_id = atom_labels.contains("alt_atom_id");

    // parse atoms
    int i_comp_id = atom_labels.at("comp_id");
    int i_atom_id = atom_labels.at("atom_id");
    int i_alt_atom_id = has_alt_atom_id ? atom_labels.at("alt_atom_id") : i_atom_id;
    int i_type_symbol = atom_labels.at("type_symbol");
    std::vector<residue::detail::Residue> residues;
    std::unordered_map<std::string, size_t> residue_names;
    for (size_t i = 0; i < atom.data.size();) { // nested loops are responsible for incrementing i
        std::string current_comp_id = atom.data[i][i_comp_id];

        residue::detail::Residue residue(current_comp_id);
        if (i+1 == atom.data.size() || atom.data[i+1][i_comp_id] != current_comp_id) {
            // single-atom residue indicates an ion
            auto element = constants::symbols::parse_element_string(atom.data[i][i_type_symbol]);
            residue.add_atom(atom.data[i][i_atom_id], constants::charge::ionic::get_charge(element), element);
            ++i;
        } else {
            // multi-atom residue
            while (i < atom.data.size() && atom.data[i][i_comp_id] == current_comp_id) {
                const std::string& atom_id = atom.data[i][i_atom_id];
                const std::string& type_symbol = atom.data[i][i_type_symbol];
                std::string alt_atom_id = atom.data[i][i_alt_atom_id];

                // sometimes the "1" in e.g. CD1 is omitted, so if the atom_id and alt_atom_id are the same anyway, add this as an alias
                if (atom_id == alt_atom_id) {
                    if (auto N = atom_id.size(); N > 2) { // second character must be a locator e.g. A, B, C, D, E, ...
                        if (atom_id[N - 1] == '1' && !std::isdigit(atom_id[N - 2])) { // last character must be 1 & previous character must be a locator
                            alt_atom_id = atom_id.substr(0, N - 1);
                        }
                    }
                }
                residue.add_atom(atom_id, alt_atom_id, constants::symbols::parse_element_string(type_symbol));
                ++i;
            }
        }

        residues.push_back(residue);
        residue_names[current_comp_id] = residues.size()-1;
    }

    // parse bonds
    if (!bond_labels.contains("comp_id") || !bond_labels.contains("atom_id_1") || !bond_labels.contains("atom_id_2") || !bond_labels.contains("value_order")) {
        throw except::io_error("CIFReader::parse_residue: Missing required labels in \"_chem_comp_bond\" section");
    }

    int i_comp_id_bond = bond_labels.at("comp_id");
    int i_atom_id_1 = bond_labels.at("atom_id_1");
    int i_atom_id_2 = bond_labels.at("atom_id_2");
    int i_value_order = bond_labels.at("value_order");
    for (size_t i = 0; i < bond.data.size();) { // nested loops are responsible for incrementing i
        const std::string& current_comp_id = bond.data[i][i_comp_id_bond];

        // parse bonds
        auto& current_residue = residues[residue_names.at(current_comp_id)];
        while (i < bond.data.size() && bond.data[i][i_comp_id_bond] == current_comp_id) {
            const std::string& atom_id_1 = bond.data[i][i_atom_id_1];
            const std::string& atom_id_2 = bond.data[i][i_atom_id_2];
            unsigned int value_order = residue::detail::Bond::parse_order(bond.data[i][i_value_order]);
            current_residue.apply_bond(residue::detail::Bond(atom_id_1, atom_id_2, value_order));
            ++i;
        }

        // remove a hydrogen bond from the N-terminus since it will almost always be bonded to the CA of the next chain
        if (const auto& m = current_residue.get_name_map(); m.contains("N")) {
            current_residue.get_atoms()[(m.at("N"))].hydrogen_bonds -= 1;
        }

        // add aliases for O and OXT that are used in e.g. GROMACS output files
        // ? it does not seem like there is any way to distinguish between the single and double bonded O's
        // ? since there is at most a single O-terminus, the impact of switching them is minimal anyway
        current_residue.add_atom("OC1", "OC1", constants::atom_t::O);
        current_residue.get_atoms().back().hydrogen_bonds = 1;
        current_residue.get_atoms().back().valency = 1; // single-bonded
        current_residue.add_atom("OC2", "OC2", constants::atom_t::O);
        current_residue.get_atoms().back().valency = 0; // double-bonded
    }

    return residues;
}

void parse_chem_comp_section(CIFSection& atom, CIFSection& bond) {
    auto residues = parse_residue(atom, bond);
    for (auto& residue : residues) {
        constants::hydrogen_atoms::residues.insert(residue.get_name(), residue.to_map());
    }
}

void parse_atom_site_section(CIFSection& atom, data::detail::AtomCollection& collection) {
    if (atom.data.empty()) {throw except::io_error("CIFReader::parse_atom_site_section: Empty data section");}
    auto labels = atom.get_label_map();

    // PDB & mmCIF equivalence:
    //   Section   	        _atom_site.group_PDB   	 
    //   Serial_No   	    _atom_site.id   	 
    //   Atom_Name   	    _atom_site.auth_atom_id   	 
    //   Alt_Loc   	        _atom_site.label_alt_id   	 
    //   Residue_Name       _atom_site.auth_comp_id   	 
    //   Strand_ID   	    _atom_site.auth_asym_id   	 
    //   Residue_No   	    _atom_site.auth_seq_id   	 
    //   Ins_Code   	    _atom_site.pdbx_PDB_ins_code   	 
    //   X   	            _atom_site.Cartn_x   	 
    //   Y   	            _atom_site.Cartn_y   	 
    //   Z   	            _atom_site.Cartn_z   	 
    //   Occupancy   	    _atom_site.occupancy   	 
    //   T_Factor   	    _atom_site.B_iso_or_equiv   	 
    //   Sigma_X   	        _atom_site.Cartn_x_esd   	 
    //   Sigma_Y   	        _atom_site.Cartn_y_esd   	 
    //   Sigma_Z   	        _atom_site.Cartn_z_esd   	 
    //   Sigma_Occupancy    _atom_site.occupancy_esd   	 
    //   Sigma_T_Factor     _atom_site.B_iso_or_equiv_esd   	 
    //   Symbol           	_atom_site.type_symbol   	 
    //   Charge           	_atom_site.pdbx_formal_charge   	 

    std::string s_residue_name = "auth_comp_id";
    {   // mandatory data 

        // prefer author labels
        if (!labels.contains(s_residue_name)) {s_residue_name = "label_comp_id";}

        if (!labels.contains("group_PDB")) {        // HETATM or ATOM
            throw except::io_error("CIFReader::parse_atom_site_section: Missing required label \"group_PDB\"");
        } if (!labels.contains(s_residue_name)) {   // residue name
            throw except::io_error("CIFReader::parse_atom_site_section: Missing required label \"comp_id\"");
        } if (!labels.contains("Cartn_x")) {        // x-coordinate
            throw except::io_error("CIFReader::parse_atom_site_section: Missing required label \"Cartn_x\"");
        } if (!labels.contains("Cartn_y")) {        // y-coordinate
            throw except::io_error("CIFReader::parse_atom_site_section: Missing required label \"Cartn_y\"");
        } if (!labels.contains("Cartn_z")) {        // z-coordinate
            throw except::io_error("CIFReader::parse_atom_site_section: Missing required label \"Cartn_z\"");
        } if (!labels.contains("type_symbol")) {    // atomic element
            throw except::io_error("CIFReader::parse_atom_site_section: Missing required label \"type_symbol\"");
        }
    }

    bool optional_data = true;
    std::string s_atom_name = "auth_atom_id";
    std::string s_chainID = "auth_asym_id";
    std::string s_residue_sequence_number = "auth_seq_id";
    {   // optional data

        // prefer author labels
        if (!labels.contains(s_atom_name)) {s_atom_name = "label_atom_id";}
        if (!labels.contains(s_residue_name)) {s_residue_name = "label_comp_id";}
        if (!labels.contains(s_chainID)) {s_chainID = "label_asym_id";}
        if (!labels.contains(s_residue_sequence_number)) {s_residue_sequence_number = "label_seq_id";}
        if (!(
            labels.contains("id") &&                        // serial number
            labels.contains("label_alt_id") &&              // alternate location
            labels.contains(s_atom_name) &&                 // atom name
            labels.contains(s_chainID) &&                   // chain ID
            labels.contains(s_residue_sequence_number) &&   // residue sequence number
            labels.contains("pdbx_PDB_ins_code") &&         // insertion code
            labels.contains("occupancy") &&                 // occupancy
            labels.contains("B_iso_or_equiv") &&            // temperature factor
            labels.contains("pdbx_formal_charge")           // charge
        )) {
            console::print_warning("CIFReader::parse_atom_site_section: Missing optional labels in \"_atom_site\" section. Non-critical data will not be loaded.");
            optional_data = false;
        }
    }

    int i_group_PDB = labels.at("group_PDB");
    int i_label_comp_id = labels.at(s_residue_name);
    int i_Cartn_x = labels.at("Cartn_x");
    int i_Cartn_y = labels.at("Cartn_y");
    int i_Cartn_z = labels.at("Cartn_z");
    int i_type_symbol = labels.at("type_symbol");

    int i_id = 0, i_label_alt_id = 0, i_label_atom_id = 0, i_label_asym_id = 0, i_label_seq_id = 0, 
        i_PDB_ins_code = 0, i_occupancy = 0, i_B_iso_or_equiv = 0, i_pdbx_formal_charge = 0;
    if (optional_data) {
        i_id = labels.at("id");
        i_label_alt_id = labels.at("label_alt_id");
        i_label_atom_id = labels.at(s_atom_name);
        i_label_asym_id = labels.at(s_chainID);
        i_label_seq_id = labels.at(s_residue_sequence_number);
        i_PDB_ins_code = labels.at("pdbx_PDB_ins_code");
        i_occupancy = labels.at("occupancy");
        i_B_iso_or_equiv = labels.at("B_iso_or_equiv");
        i_pdbx_formal_charge = labels.at("pdbx_formal_charge");
    }

    auto shorten = [] (const std::string& s) -> std::string {
        return s.size() < 6 ? s : s.substr(0, 7);
    };

    int discarded_hydrogens = 0;
    for (size_t i = 0; i < atom.data.size(); ++i) {
        auto& group_PDB = atom.data[i][i_group_PDB];
        if (data::record::Record::get_type(group_PDB) != data::record::RecordType::ATOM) {
            throw except::io_error("CIFReader::parse_atom_site_section: Unrecognized group_PDB \"" + group_PDB + "\"");
        }

        int serial = 0, resSeq = 0;
        double occupancy = 1, tempFactor = 0;
        char chainID = ' ';
        std::string name, altLoc, resName, iCode, charge;
        Vector3<double> coords;
        constants::atom_t element;

        // load mandatory data
        try {
            name = atom.data[i][i_label_atom_id];
            resName = atom.data[i][i_label_comp_id];
            coords = {
                std::stod(shorten(atom.data[i][i_Cartn_x])), 
                std::stod(shorten(atom.data[i][i_Cartn_y])), 
                std::stod(shorten(atom.data[i][i_Cartn_z])),
            };
            element = constants::symbols::parse_element_string(atom.data[i][i_type_symbol]);
            if (optional_data) {
                altLoc = atom.data[i][i_label_alt_id];
                chainID = atom.data[i][i_label_asym_id][0];
                iCode = atom.data[i][i_PDB_ins_code];
                charge = atom.data[i][i_pdbx_formal_charge];
                if (!atom.data[i][i_id].starts_with('.')) {serial = std::stoi(atom.data[i][i_id]);}
                if (!atom.data[i][i_label_seq_id].starts_with('.')) {resSeq = std::stoi(atom.data[i][i_label_seq_id]);}
                if (!atom.data[i][i_occupancy].starts_with('.')) {occupancy = std::stod(shorten(atom.data[i][i_occupancy]));}
                if (!atom.data[i][i_B_iso_or_equiv].starts_with('.')) {tempFactor = std::stod(shorten(atom.data[i][i_B_iso_or_equiv]));}
            }
        } catch (const std::exception& e) {
            console::print_warning(
                "CIFReader::parse_atom_site_section: Invalid field values in line: \n\"" + 
                std::accumulate(
                    atom.data[i].begin(), atom.data[i].end(), std::string(), 
                    [] (const std::string& a, const std::string& b) {return a + " " + b;}
                ) + "\".");
            throw e;
        }
        data::record::Atom a(serial, name, altLoc, resName, chainID, resSeq, iCode, coords, occupancy, tempFactor, element, charge);

        // check if this is a hydrogen atom
        if (a.element == constants::atom_t::H && !settings::general::keep_hydrogens) {
            discarded_hydrogens++;
            continue;
        }

        // check if this is a water molecule
        if (a.is_water()) {collection.add(data::record::Water(std::move(a)));} 
        else {collection.add(std::move(a));}
    }

    if (!settings::molecule::use_occupancy) {
        for (auto& a : collection.protein_atoms) {a.occupancy = 1.0;}
    }

    if (discarded_hydrogens != 0) {
        console::print_text("Discarded " + std::to_string(discarded_hydrogens) + " explicit hydrogen atoms.");
    }
}

// Extract the current section from the file. This also advances the input stream to the next section
CIFSection extract_section(std::string line, std::ifstream& input) {
    std::vector<std::string> labels;
    std::vector<std::vector<std::string>> data;

    // read labels
    data.push_back({}); // add an empty data entry in case the labels also contains data
    do {
        line = utility::remove_all(line, "\n\r");
        if (line.starts_with('#')) {continue;}
        if (line.starts_with('_')) {
            auto tokens = utility::split(line, ' ');
            labels.push_back(utility::split(tokens[0], '.').back());

            // if the label is followed by a value, add it to the data
            if (tokens.size() == 2) {data[0].push_back(tokens[1]);}
            continue;
        }
    } while(input.peek() == '_' && getline(input, line));

    // read data
    if (data[0].empty()) {data.pop_back();} // remove the dummy data if empty
    while(input.peek() != '_' && getline(input, line) && !line.starts_with("loop_")) {
        if (line.starts_with('#')) {continue;}
        auto values = utility::split(line, " ");

        int concatenated = 0;
        for (size_t i = 0; i < values.size(); i++) {
            if (values[i].starts_with('"') && !values[i].ends_with('"')) {
                if (i < values.size() && values[i+1].ends_with('"')) {
                    values[i] = values[i].substr(1) + " " + values[i+1].substr(0, values[i+1].size()-1);
                } else {
                    throw except::io_error("CIFReader::read: Unterminated quote in data section: \n\"" + line + "\"");
                }
                ++concatenated;
                continue;
            }
            values[i-concatenated] = values[i];
        }
        values.resize(values.size()-concatenated);

        // sanity check
        if (values.size() != labels.size()) {
            throw except::io_error(
                "CIFReader::read: Inconsistent number of values in data section: \n\"" + line + "\"" +
                "\nExpected " + std::to_string(labels.size()) + " values, got " + std::to_string(values.size())
            );
        }

        data.push_back(values);
    };

    return CIFSection{labels, data};
}

std::vector<residue::detail::Residue> io::detail::CIFReader::read_residue(const io::File& path) {
    std::ifstream input(path);
    if (!input.is_open()) {throw except::io_error("CIFReader::read_residue: Could not open file \"" + path.str() + "\"");}

    std::string line;
    CIFSection chem_comp_atom, chem_comp_bond;
    while(getline(input, line)) {
        if (line.find("_chem_comp_atom.") != std::string::npos) {
            chem_comp_atom = extract_section(line, input);
        }

        if (line.find("_chem_comp_bond.") != std::string::npos) {
            chem_comp_bond = extract_section(line, input);
        }
    }

    if (chem_comp_bond.empty()) {
        if (chem_comp_atom.empty()) {throw except::io_error("CIFReader::read_residue: Could not find any data sections in file \"" + path.str() + "\"");}
        return {parse_ion(chem_comp_atom)};
    } else {
        return parse_residue(chem_comp_atom, chem_comp_bond);
    }
}

void io::detail::CIFReader::read(const io::File& path) {
    console::print_info("Reading CIF file from \"" + path.str() + "\"");
    console::indent();

    std::ifstream input(path);
    if (!input.is_open()) {throw except::io_error("CIFReader::read: Could not open file \"" + path.str() + "\"");}

    std::string line;
    CIFSection chem_comp_atom, chem_comp_bond, atom_site;
    while(getline(input, line)) {
        if (line.find("_chem_comp_atom.") != std::string::npos) {
            chem_comp_atom = extract_section(line, input);
        }

        if (line.find("_chem_comp_bond.") != std::string::npos) {
            chem_comp_bond = extract_section(line, input);
        }

        if (line.find("_atom_site.") != std::string::npos) {
            atom_site = extract_section(line, input);
        }
    }

    if (atom_site.empty()) {throw except::io_error("CIFReader::read: Could not find any atomic data section in file \"" + path.str() + "\"");}
    if (!chem_comp_atom.empty() && !chem_comp_bond.empty()) {parse_chem_comp_section(chem_comp_atom, chem_comp_bond);}
    parse_atom_site_section(atom_site, *file);

    unsigned int n_pa = file->protein_atoms.size();
    unsigned int n_ha = file->hydration_atoms.size();

    console::print_text("Successfully read " + std::to_string(n_pa + n_ha) + " atomic records.");
    if (n_ha != 0) {console::print_text("\t" + std::to_string(file->hydration_atoms.size()) + " of these are hydration atoms.");}
    console::unindent();
}