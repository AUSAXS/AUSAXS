// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <residue/ResidueStorage.h>

#include <constants/Constants.h>
#include <io/ExistingFile.h>
#include <io/File.h>
#include <residue/detail/InvalidResidueMap.h>
#include <residue/detail/ResidueMap.h>
#include <residue/detail/ResidueStorageBasis.h>
#include <settings/GeneralSettings.h>
#include <settings/MoleculeSettings.h>
#include <utility/Console.h>
#include <utility/Curl.h>
#include <utility/StringUtils.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::residue;
using namespace ausaxs::residue::detail;

namespace {
    struct AtomAlias {
        std::string residue; // empty: any residue
        std::string name;
    };

    // The names other conventions give an atom, with the CCD name the residue definitions carry it under.
    // The two C-terminal oxygens are interchangeable; which of them carries the hydrogen is a convention.
    const std::unordered_multimap<std::string, AtomAlias> atom_aliases = {
        {"oc1", {.residue="",    .name="oxt"}}, // GROMACS
        {"oc2", {.residue="",    .name="o"}},
        {"ot1", {.residue="",    .name="oxt"}}, // CHARMM
        {"ot2", {.residue="",    .name="o"}},
        {"cd",  {.residue="ILE", .name="cd1"}}
    };

    /**
     * @brief The name the residue definition carries this atom under. A name the definition has itself is never replaced by an alias.
     */
    std::string resolve_atom_name(const ResidueMap& residue, const std::string& residue_name, const std::string& atom_name, constants::atom_t atom) {
        if (residue.contains(atom_name, atom)) {return atom_name;}
        auto [first, last] = atom_aliases.equal_range(utility::to_lowercase(atom_name));
        for (; first != last; ++first) {
            const auto& alias = first->second;
            if (!alias.residue.empty() && alias.residue != residue_name) {continue;}
            if (residue.contains(alias.name, atom)) {return alias.name;}
        }
        return atom_name;
    }

    struct ProtonationSite {
        std::string atom;
        constants::atom_t element;
        int hydrogens;
    };

    struct ProtonationVariant {
        std::string name, base;
        std::vector<ProtonationSite> sites; // absolute counts, so a variant does not depend on the state of its standard residue
    };

    // AMBER and CHARMM name a residue by its protonation state. Most of these names are also CCD codes for unrelated compounds 
    // (HSE is homoserine, HIP phosphohistidine), so they are only read as such when settings::molecule::force_field_residue_names is set.
    const std::vector<ProtonationVariant>& protonation_variant_table() {
        constexpr auto N = constants::atom_t::N, O = constants::atom_t::O, S = constants::atom_t::S;
        static const std::vector<ProtonationVariant> variants = {
            {"HID", "HIS", {{"nd1", N, 1}, {"ne2", N, 0}}}, {"HSD", "HIS", {{"nd1", N, 1}, {"ne2", N, 0}}}, // hydrogen on ND1 only
            {"HIE", "HIS", {{"nd1", N, 0}, {"ne2", N, 1}}}, {"HSE", "HIS", {{"nd1", N, 0}, {"ne2", N, 1}}}, // hydrogen on NE2 only
            {"HIP", "HIS", {{"nd1", N, 1}, {"ne2", N, 1}}}, {"HSP", "HIS", {{"nd1", N, 1}, {"ne2", N, 1}}}, // both
            {"ASH", "ASP", {{"od1", O, 0}, {"od2", O, 1}}}, {"GLH", "GLU", {{"oe1", O, 0}, {"oe2", O, 1}}}, // neutral acids
            {"CYX", "CYS", {{"sg",  S, 0}}},                {"CYM", "CYS", {{"sg",  S, 0}}},                // disulfide-bonded and deprotonated cysteine
            {"LYN", "LYS", {{"nz",  N, 2}}}                                                                 // neutral lysine
        };
        return variants;
    }

    void warn_force_field_name(const std::string& name) {
        static std::unordered_set<std::string> warned;
        if (!warned.insert(name).second) {return;}
        auto it = std::ranges::find(protonation_variant_table(), name, &ProtonationVariant::name);
        console::print_warning(
            "Residue \"" + name + "\" is read by its CCD definition, but AMBER and CHARMM also use this name for a protonation state of " + it->base + ". "
            "If the structure uses force-field names, read them as such with 'hydrogens --force-field-names', "
            "or keep the explicit hydrogens of the file with 'hydrogens --keep' for higher precision."
        );
    }
}

ResidueStorage::ResidueStorage() = default;

ResidueStorage::~ResidueStorage() = default;

void ResidueStorage::insert(const std::string& name, ResidueMap&& residue) {
    data.emplace(name, std::move(residue));
}

bool ResidueStorage::contains(const std::string& name) {
    if (!initialized) {initialize();}
    if (settings::molecule::force_field_residue_names && protonation_variants.contains(name)) {return true;}
    return data.contains(name);
}

bool ResidueStorage::load(const std::string& name) {
    if (!initialized) {initialize();}
    if (protonation_variants.contains(name)) {
        if (settings::molecule::force_field_residue_names) {return true;}
        warn_force_field_name(name);
    }
    if (data.contains(name)) {return true;}

    // small cache to avoid spamming the console with the same download
    static std::unordered_map<std::string, bool> seen_before;
    if (seen_before.contains(name)) {return false;}
    seen_before[name] = true;

    console::indent(2);
    console::print_text("Unknown residue: \"" + name + "\". Attempting to download specification.");
    console::indent();
    bool downloaded = false;
    try {
        downloaded = update_or_download_residue(name);
    } catch (const except::io_error& e) {
        console::print_warning(e.what());
    }
    console::unindent(3);
    return downloaded;
}

ResidueMap& ResidueStorage::get(const std::string& name) {
    if (!load(name)) {
        // singleton for all unknown residues which will return 0 hydrogens for all atoms
        static InvalidResidueMap invalid_residue;
        return invalid_residue;
    }
    if (settings::molecule::force_field_residue_names) {
        if (auto it = protonation_variants.find(name); it != protonation_variants.end()) {return it->second;}
    }
    return data.at(name);
}

void ResidueStorage::initialize() {
    initialized = true;
    const std::string path = settings::general::helper::residue_folder();
    io::Folder(path).create();
    std::ifstream file(path + "master.dat");
    if (!file.is_open()) {
        write_master_basis();
        file.open(path + "master.dat");
        if (!file.is_open()) {throw except::io_error("ResidueStorage::initialize: Could not open file: " + path + "master" + ".dat");}
    }

    std::string line;
    while (file.peek() != EOF) {
        std::getline(file, line);

        // skip until we reach the start of a residue
        if (line.find('#') == std::string::npos) {
            continue;
        }
        // the line following the # is the name of the residue
        std::getline(file, line);
        std::string residue = line;

        // prepare map
        std::unordered_map<AtomKey, int> map;
        while (file.peek() != EOF) {
            std::getline(file, line);
            // stop if we reach the start of a new residue
            if (line.empty() || line.find('#') != std::string::npos) {
                break;
            }

            // lines are of the form "element atom hydrogens"
            std::vector<std::string> tokens = utility::split(line, " \n\r");
            if (tokens.size() != 3) {throw except::io_error("ResidueStorage::initialize: Invalid line in master file: " + line + ". Perhaps the file is corrupted. Delete it to regenerate.");}
            const std::string& element = tokens[0];
            const std::string& atom = tokens[1];
            int hydrogens = std::stoi(tokens[2]);
            map.emplace(AtomKey(atom, constants::symbols::parse_element_string(element)), hydrogens);
        }
        insert(residue, std::move(map));
    }
    insert_protonation_variants();
}

void ResidueStorage::insert_protonation_variants() {
    for (const auto& variant : protonation_variant_table()) {
        if (!data.contains(variant.base)) {continue;}
        auto map = data.at(variant.base).get_backing_map();
        for (const auto& site : variant.sites) {
            map.insert_or_assign(AtomKey(site.atom, site.element), site.hydrogens);
        }
        protonation_variants.insert_or_assign(variant.name, ResidueMap(std::move(map)));
    }
}

bool ResidueStorage::update_or_download_residue(const std::string& name) {
    const std::string path = settings::general::helper::residue_folder();
    std::regex regex("[A-Z0-9]{1,3}");

    if (std::regex_match(name, regex)) {
        // check if the file already exists. if not, download it.
        if (!std::filesystem::exists(path + name + ".cif")) {
            if (settings::general::offline) {
                console::print_warning("ResidueStorage::download_residue: \"" + name + "\" specification cannot be downloaded as offline mode is enabled!");
                return false;
            }

            bool downloaded = curl::download("files.rcsb.org/ligands/view/" + name + ".cif", path + name + ".cif"); // download the cif file
            if (!downloaded) {return false;}
        } else {
            console::print_text_minor("Residue " + name + " is already downloaded, but not present in the master list.");
            console::print_text_minor("The file will now be parsed and re-added to the master file.");
        }

        // parse the cif file & add it to storage
        residue::detail::ResidueMap map;
        try {
            map = Residue::parse(path + name + ".cif").to_map();
        } catch (const std::exception& e) {
            // if the residue could not be parsed, try to redownload it and parse it again
            console::print_warning("ResidueStorage::download_residue: Could not parse residue: " + name + ". \nError: " + e.what());
            console::print_text("Attempting to redownload the residue definition.");

            // the cached file is unusable; remove it so a failed redownload does not leave it to fail again on every later run
            io::File(path + name + ".cif").remove();
            if (!curl::download("files.rcsb.org/ligands/view/" + name + ".cif", path + name + ".cif")) {return false;}
            map = Residue::parse(path + name + ".cif").to_map();
        }
        insert(name, std::move(map));

        // write the residue to the master file
        write_residue(name);
    } else {
        throw except::io_error("ResidueStorage::download_residue: Invalid residue name: \"" + name + "\"; expected a 2-3 letter code.");
    }
    return true;
}

void ResidueStorage::write_residue(const std::string& name) {
    const std::string path = settings::general::helper::residue_folder();
    io::Folder(path).create();
    std::ofstream file(path + "master.dat", std::ios::app); // open in append mode
    if (!file.is_open()) {throw except::io_error("ResidueStorage::write_residue: Could not open file: " + path + "master" + ".dat");}

    // write the map to the master file
    auto map = get(name);
    file << "#" << "\n" << name << "\n"; // residue header
    for (const auto& [key, val] : map.get_backing_map()) {
        file << constants::symbols::to_string(key.atom) << " " << key.name << " " << val << "\n";
    }
    file << std::endl;
}

constants::atomic_group_t ResidueStorage::get_atomic_group(const std::string& residue_name, const std::string& atom_name, constants::atom_t atom) {
    auto& residue = get(residue_name);
    return residue.get_atomic_group(resolve_atom_name(residue, residue_name, atom_name, atom), atom);
}

double ResidueStorage::get_hydrogens(const std::string& residue_name, const std::string& atom_name, constants::atom_t atom) {
    auto& residue = get(residue_name);
    return residue.get(resolve_atom_name(residue, residue_name, atom_name, atom), atom);
}