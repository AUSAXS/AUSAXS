#pragma once

#include <map>
#include <string>
#include <iostream>

#include <utility/Exceptions.h>
#include <utility/SimpleResidueMap.h>

saxs::detail::SimpleResidueMap::SimpleResidueMap() : SimpleMap() {}

saxs::detail::SimpleResidueMap::SimpleResidueMap(std::unordered_map<std::string, unsigned int> names, std::unordered_map<std::string, std::string> elements) : SimpleMap(names) {
    calculate_average(elements);
}

const unsigned int& saxs::detail::SimpleResidueMap::get(AtomKey key) const {
    // first check if the key is in the map
    if (data.find(key.name) != data.end()) {return data.at(key.name);}

    // if not, check if the key is a hydrogen
    if (key.symbol == "H") {return 0;}

    // estimate the number of bonds as the average for that element
    if (average.find(key.symbol) != average.end()) {
        return average.at(key.symbol);
    } else {
        throw except::map_error("Error in SimpleResidueMap::get: Key " + key.name + " not found in map, and no estimate for element " + key.symbol + " is available.");
    }
}

void saxs::detail::SimpleResidueMap::calculate_average(std::unordered_map<std::string, std::string> elements) {
    std::unordered_map<std::string, unsigned int> counts;
    
    for (auto& [key, value] : data) {
        if (elements.find(key) == elements.end()) {
            throw except::map_error("Error in SimpleResidueMap::calculate_average: Key " + key + " not found in element map");
        }
        std::string symbol = elements.at(key);
        average[symbol] += value;
        counts[symbol] += 1;
    }

    for (auto& [key, value] : average) {
        average[key] /= counts[key];
    }
}
