/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <residue/ResidueMap.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <constants/Constants.h>
#include <settings/MoleculeSettings.h>

#include <string>

using namespace residue::detail;

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
        throw except::map_error("SimpleResidueMap::get: Key " + key.name + " not found in map, and no estimate for element id " + constants::symbols::to_string(key.atom) + " is available.");
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

const std::unordered_map<AtomKey, int>& ResidueMap::get_map() const {return map;}

constants::atomic_group_t ResidueMap::get_atomic_group(const std::string& atom_name, constants::atom_t atom_type) {
    auto key = AtomKey(atom_name, atom_type);
    if (!map.contains(key)) {
        if (key.atom == constants::atom_t::H) {return constants::atomic_group_t::unknown;}
        if (settings::molecule::throw_on_unknown_atom) {
            throw except::map_error("ResidueMap::get_atomic_group: Key " + atom_name + " not found in map.");
        } else {
            return constants::atomic_group_t::unknown;
        }
    }
    int hydrogens = map.at(key);
    return constants::symbols::get_atomic_group(atom_type, hydrogens);
}
