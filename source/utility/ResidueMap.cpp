#include <string>
#include <iostream>

#include <utility/Exceptions.h>
#include <utility/ResidueMap.h>

using namespace saxs::detail;

ResidueMap::ResidueMap() {}

ResidueMap::ResidueMap(std::unordered_map<AtomKey, unsigned int> map) {
    this->map = map;
    this->calculate_average();
}

double ResidueMap::get(AtomKey key) {
    // first check if the key is in the map
    if (map.find(key) != map.end()) {return map.at(key);}

    // if not, check if the key is a hydrogen
    if (key.symbol == "H") {return 0;}

    // estimate the number of bonds as the average for that element
    if (__builtin_expect(update_average, false)) {this->calculate_average();}
    if (average.find(key.symbol) != average.end()) {
        return average.at(key.symbol);
    } else {
        throw except::map_error("Error in SimpleResidueMap::get: Key " + key.name + " not found in map, and no estimate for element " + key.symbol + " is available.");
    }
}

void ResidueMap::insert(AtomKey key, unsigned int value) {
    map[key] = value;
    update_average = true;
}

void ResidueMap::insert(std::string name, std::string symbol, unsigned int value) {
    insert(AtomKey(name, symbol), value);
}

void ResidueMap::calculate_average() {
    std::unordered_map<std::string, unsigned int> counts;
    
    for (auto& [key, value] : map) {
        average[key.symbol] += value;
        counts[key.symbol] ++;
    }

    for (auto& [key, value] : average) {
        average[key] /= counts[key];
    }
}
