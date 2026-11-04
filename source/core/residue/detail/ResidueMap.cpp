// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <residue/detail/ResidueMap.h>

#include <constants/Constants.h>
#include <settings/MoleculeSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

#include <string>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::residue::detail;

AtomKey::AtomKey(const std::string& name, constants::atom_t atom) : name(utility::to_lowercase(name)), atom(atom) {}
bool AtomKey::operator==(const AtomKey& other) const {
    return name == other.name;
}

std::size_t std::hash<residue::detail::AtomKey>::operator()(const AtomKey& k) const {return std::hash<std::string>()(k.name);}

ResidueMap::ResidueMap() = default;

ResidueMap::ResidueMap(std::unordered_map<AtomKey, int> map) : map(std::move(map)) {
    this->calculate_average();
}

namespace {
    /**
     * @brief The name a force field gives an atom that the residue definitions know by another name, or an empty string if there is none.
     *        CHARMM calls the two C-terminal oxygens OT1 and OT2, where the definitions carry the GROMACS names OC1 and OC2.
     */
    std::string alias_of(const std::string& lowercase_name) {
        if (lowercase_name == "ot1") {return "oc1";}
        if (lowercase_name == "ot2") {return "oc2";}
        return {};
    }
}

const int* ResidueMap::find(const AtomKey& key) const {
    if (auto it = map.find(key); it != map.end()) {return &it->second;}
    if (auto alias = alias_of(key.name); !alias.empty()) {
        if (auto it = map.find(AtomKey(alias, key.atom)); it != map.end()) {return &it->second;}
    }
    return nullptr;
}

bool ResidueMap::contains(const std::string& name, constants::atom_t atom) const {
    return map.contains(AtomKey(name, atom));
}

double ResidueMap::get(const AtomKey& key) {
    // first check if the key is in the map
    if (auto value = find(key)) {return *value;}

    // if not, check if the key is a hydrogen
    if (key.atom == constants::atom_t::H) {return 0;}

    // estimate the number of bonds as the average for that element
    if (update_average) [[unlikely]] {this->calculate_average();}
    if (average.contains(key.atom)) {
        return average.at(key.atom);
    }
    
    if (!settings::molecule::allow_unknown_atoms) {
        throw except::map_error("ResidueMap::get: Key " + key.name + " not found in map, and no estimate for element id " + constants::symbols::to_string(key.atom) + " is available.");
    }

    static bool warned = false;
    if (!warned) {
        console::print_warning(
            "ResidueMap::get: Key " + key.name + " not found in map, and no estimate for element id " + constants::symbols::to_string(key.atom) + " is available."
            "Further warnings of this type will be suppressed."
        );
        warned = true;
    }
    return 0;
}

void ResidueMap::insert(const AtomKey& key, int value) {
    map[key] = value;
    update_average = true;
}

void ResidueMap::insert(const std::string& name, constants::atom_t symbol, int value) {
    insert(AtomKey(name, symbol), value);
}

void ResidueMap::calculate_average() {
    std::unordered_map<constants::atom_t, int> counts;
    
    for (auto& [key, value] : map) {
        average[key.atom] += value;
        counts[key.atom] ++;
    }

    for (auto& [key, value] : average) {
        average[key] /= counts[key];
    }
}

std::unordered_map<AtomKey, int>& ResidueMap::get_backing_map() {return map;}

constants::atomic_group_t ResidueMap::get_atomic_group(const std::string& atom_name, constants::atom_t atom_type) {
    auto key = AtomKey(atom_name, atom_type);
    auto value = find(key);
    if (!value) {
        if (key.atom == constants::atom_t::H) {return constants::atomic_group_t::unknown;}
        if (!settings::molecule::allow_unknown_atoms) {
            throw except::map_error("ResidueMap::get_atomic_group: Key " + atom_name + " not found in map.");
        }
        static bool warned = false;
        if (!warned) {
            console::print_warning(
                "ResidueMap::get_atomic_group: Key " + atom_name + " not found in map."
                "Further warnings of this type will be suppressed."
            );
            warned = true;
        }
        return constants::atomic_group_t::unknown;
    }
    int hydrogens = *value;
    return constants::symbols::get_atomic_group(atom_type, hydrogens);
}
