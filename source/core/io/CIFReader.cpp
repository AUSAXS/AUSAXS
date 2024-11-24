#include <io/CIFReader.h>
#include <io/File.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <residue/ResidueParser.h>
#include <constants/Constants.h>

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

void parse_ion(const CIFSection& atom) {
    if (atom.data.size() != 1) {throw except::io_error("CIFReader::parse_ion: Invalid number of data entries in atom section");}
    auto labels = atom.get_label_map();

    if (!labels.contains("formula") || !labels.contains("pdbx_formal_charge")) {
        throw except::io_error("CIFReader::parse_ion: Missing required labels in \"_chem_comp\" section");
    }

    std::string formula = atom.data[0][labels.at("formula")];
    if (constants::hydrogen_atoms::residues.contains(formula)) {return;}

    if (!constants::symbols::detail::string_to_atomt_map.contains(formula)) {
        throw except::io_error("CIFReader::parse_ion: Unrecognized ion: \"" + formula + "\"");
    }

    int charge = std::stoi(atom.data[0][labels.at("pdbx_formal_charge")]);
    residue::detail::Residue residue(formula);
    residue.add_atom(formula, charge, constants::symbols::parse_element_string(formula));
    constants::hydrogen_atoms::residues.insert_and_write(formula, residue.to_map());
}

void parse_residue(CIFSection& atom, CIFSection& bond) {
    if (atom.data.empty() || bond.data.empty()) {throw except::io_error("CIFReader::parse_residue: Empty data section");}
    auto atom_labels = atom.get_label_map();
    auto bond_labels = bond.get_label_map();

    if (!atom_labels.contains("comp_id") || !atom_labels.contains("atom_id") || !atom_labels.contains("alt_atom_id")) {
        throw except::io_error("CIFReader::parse_residue: Missing required labels in \"_chem_comp_atom\" section");
    }

    // parse atoms
    int i_comp_id = atom_labels.at("comp_id");
    int i_atom_id = atom_labels.at("atom_id");
    int i_alt_atom_id = atom_labels.at("alt_atom_id");
    int i_type_symbol = atom_labels.at("type_symbol");
    std::vector<residue::detail::Residue> residues;
    std::unordered_map<std::string, size_t> residue_names;
    for (size_t i = 0; i < atom.data.size(); ++i) {
        std::string current_comp_id = atom.data[i][i_comp_id];

        // if the residue is already present, skip it
        if (constants::hydrogen_atoms::residues.contains(current_comp_id)) {
            while (i < atom.data.size() && atom.data[i][i_comp_id] == current_comp_id) {++i;}
        }

        // else, parse all atoms in the residue
        residue::detail::Residue residue(current_comp_id);
        while (i < atom.data.size() && atom.data[i][i_comp_id] == current_comp_id) {
            auto& atom_id = atom.data[i][i_atom_id];
            auto& alt_atom_id = atom.data[i][i_alt_atom_id];
            auto& type_symbol = atom.data[i][i_type_symbol];

            // Sometimes the "1" in e.g. CD1 is omitted
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
        if (residue.get_atoms().size() == 1) {throw except::io_error("CIFReader::parse_residue: Residue \"" + current_comp_id + "\" has only one atom (ion?)");}
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
    for (size_t i = 0; i < bond.data.size(); ++i) {
        std::string current_comp_id = bond.data[i][i_comp_id_bond];

        // skip if not present in atom section
        if (!residue_names.contains(current_comp_id)) {
            while (i < bond.data.size() && bond.data[i][i_comp_id_bond] == current_comp_id) {++i;}
        }

        // parse bonds
        auto& current_residue = residues[residue_names.at(current_comp_id)];
        while (i < bond.data.size() && bond.data[i][i_comp_id_bond] == current_comp_id) {
            auto& atom_id_1 = bond.data[i][i_atom_id_1];
            auto& atom_id_2 = bond.data[i][i_atom_id_2];
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

    for (auto& residue : residues) {
        constants::hydrogen_atoms::residues.insert_and_write(residue.get_name(), residue.to_map());
    }
}

void parse_chem_comp_section(CIFSection& atom, CIFSection& bond) {
    if (bond.empty()) {parse_ion(atom);}
    else {parse_residue(atom, bond);}
}

void parse_atom_site_section(CIFSection&) {

}

void io::detail::CIFReader::read(const io::File& path) {
    // check if file was succesfully opened
    std::ifstream input(path);
    if (!input.is_open()) {throw except::io_error("CIFReader::read: Could not open file \"" + path.str() + "\"");}

    // Extract the current section from the file. This also advances the input stream to the next section
    auto extract_section = [&] () -> CIFSection {
        std::vector<std::string> labels;
        std::vector<std::vector<std::string>> data;
        std::string line;

        // read labels
        data.push_back({}); // add an empty data entry in case the labels also contains data
        while(getline(input, line) && line.find("loop_") == std::string::npos) {
            line = utility::remove_all(line, "\n\r");
            if (line.starts_with('#')) {continue;}
            if (line.starts_with('_')) {
                auto tokens = utility::split(line, ' ');
                labels.push_back(utility::split(line, '.').back());

                // if the label is followed by a value, add it to the data
                if (tokens.size() == 2) {
                    data[0].push_back(tokens[1]);
                }
                continue;
            }
            break;
        }

        // read data
        if (data[0].empty()) {data.pop_back();} // remove the data if empty
        while(getline(input, line) && line.find("loop_") == std::string::npos) {
            auto values = utility::split(line, " ");
            if (values.size() != labels.size()) {
                throw except::io_error("CIFReader::read: Inconsistent number of values in data section: \n\"" + line + "\"");
            }

            int concatenated = 0;
            for (size_t i = 0; i < values.size(); i++) {
                if (values[i].starts_with('"')) {
                    if (i < values.size() && values[i+1].ends_with('"')) {
                        values[i] = values[i].substr(1) + " " + values[i+1].substr(0, values[i+1].size()-1);
                        continue;
                    } else {
                        throw except::io_error("CIFReader::read: Unterminated quote in data section: \n\"" + line + "\"");
                    }
                    ++concatenated;
                    continue;
                }
                values[i-concatenated] = values[i];
            }
            values.resize(values.size()-concatenated);
            data.push_back(values);
        }
        return CIFSection{labels, data};
    };

    std::string line;
    CIFSection chem_comp_atom, chem_comp_bond, atom_site;
    while(getline(input, line)) {
        if (line.find("_chem_comp_atom") != std::string::npos) {
            chem_comp_atom = extract_section();
        }

        if (line.find("_chem_comp_bond") != std::string::npos) {
            chem_comp_bond = extract_section();
        }

        if (line.find("_atom_site") != std::string::npos) {
            chem_comp_atom = extract_section();
        }
    }

    if (atom_site.empty() && chem_comp_atom.empty() && chem_comp_bond.empty()) {throw except::io_error("CIFReader::read: Could not find any data sections in file \"" + path.str() + "\"");}
    if (!atom_site.empty()) {parse_atom_site_section(atom_site);}
    if (int v = chem_comp_atom.empty() + chem_comp_bond.empty(); v == 1) {throw except::io_error("CIFReader::read: Incomplete data sections in file \"" + path.str() + "\"");}
    else if (v != 0) {parse_chem_comp_section(chem_comp_atom, chem_comp_bond);}
}