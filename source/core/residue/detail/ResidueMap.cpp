// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <residue/detail/ResidueMap.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <constants/Constants.h>
#include <settings/MoleculeSettings.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::residue::detail;

AtomKey::AtomKey(const std::string& name, constants::atom_t atom) : name(utility::to_lowercase(name)), atom(atom) {}
bool AtomKey::operator==(const AtomKey& other) const {
    return name == other.name;
}

unsigned int std::hash<residue::detail::AtomKey>::operator()(const AtomKey& k) const {return std::hash<std::string>()(k.name);}

ResidueMap::ResidueMap() = default;

ResidueMap::ResidueMap(const std::unordered_map<AtomKey, int>& map) {
    this->map = map;
    this->calculate_average();
}

bool ResidueMap::contains(const std::string& name, constants::atom_t atom) const {
    return map.contains(AtomKey(name, atom));
}

double ResidueMap::get(const AtomKey& key) {
    // first check if the key is in the map
    if (map.contains(key)) {return map.at(key);}

    // if not, check if the key is a hydrogen
    if (key.atom == constants::atom_t::H) {return 0;}

    // estimate the number of bonds as the average for that element
    if (update_average) [[unlikely]] {this->calculate_average();}
    if (average.contains(key.atom)) {
        return average.at(key.atom);
    } else {
        if (settings::molecule::throw_on_unknown_atom) {
            throw except::map_error("ResidueMap::get: Key " + key.name + " not found in map, and no estimate for element id " + constants::symbols::to_string(key.atom) + " is available.");
        } else {
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
    }
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
    if (!map.contains(key)) {
        if (key.atom == constants::atom_t::H) {return constants::atomic_group_t::unknown;}
        if (settings::molecule::throw_on_unknown_atom) {
            throw except::map_error("ResidueMap::get_atomic_group: Key " + atom_name + " not found in map.");
        } else {
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
    }
    int hydrogens = map.at(key);
    return constants::symbols::get_atomic_group(atom_type, hydrogens);
}
