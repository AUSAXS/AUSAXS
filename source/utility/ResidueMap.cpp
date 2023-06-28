#include <string>
#include <iostream>

#include <utility/Exceptions.h>
#include <utility/ResidueMap.h>
#include <utility/StringUtils.h>

using namespace saxs::detail;

saxs::detail::AtomKey::AtomKey(const std::string& name, const std::string& symbol) : name(utility::to_lowercase(name)), symbol(symbol) {}
bool saxs::detail::AtomKey::operator==(const AtomKey& other) const {
    return name == other.name;
}

unsigned int std::hash<saxs::detail::AtomKey>::operator()(const saxs::detail::AtomKey& k) const {return std::hash<std::string>()(k.name);}

ResidueMap::ResidueMap() = default;

ResidueMap::ResidueMap(const std::unordered_map<AtomKey, int>& map) {
    this->map = map;
    this->calculate_average();
}

double ResidueMap::get(const AtomKey& key) {
    // first check if the key is in the map
    if (map.contains(key)) {return map.at(key);}

    // if not, check if the key is a hydrogen
    if (key.symbol == "H") {return 0;}

    // estimate the number of bonds as the average for that element
    if (update_average) [[unlikely]] {this->calculate_average();}
    if (average.contains(key.symbol)) {
        return average.at(key.symbol);
    } else {
        throw except::map_error("SimpleResidueMap::get: Key " + key.name + " not found in map, and no estimate for element " + key.symbol + " is available.");
    }
}

void ResidueMap::insert(const AtomKey& key, int value) {
    map[key] = value;
    update_average = true;
}

void ResidueMap::insert(const std::string& name, const std::string& symbol, int value) {
    insert(AtomKey(name, symbol), value);
}

void ResidueMap::calculate_average() {
    std::unordered_map<std::string, int> counts;
    
    for (auto& [key, value] : map) {
        average[key.symbol] += value;
        counts[key.symbol] ++;
    }

    for (auto& [key, value] : average) {
        average[key] /= counts[key];
    }
}

const std::unordered_map<AtomKey, int>& ResidueMap::get_map() const {return map;}