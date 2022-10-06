#pragma once

#include <utility/Utility.h>

#include <unordered_map>
#include <string>

namespace saxs {
    namespace detail {
        class SimpleResidueMap {
            public:
                /**
                 * @brief The key type used in the SimpleResidueMap. 
                 *        The symbol is required to avoid ambiguities since the name is always capitalized in PDB files,
                 *        so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca). 
                 */
                struct AtomKey {
                    AtomKey(std::string name, std::string symbol) : name(utility::to_lowercase(name)), symbol(symbol) {}
                    std::string name;
                    std::string symbol;
                };

                /**
                 * @brief Create a new empty SimpleResidueMap.
                 */
                SimpleResidueMap();

                /**
                 * @brief Create a new SimpleResidueMap from a std::map.
                 */
                SimpleResidueMap(std::unordered_map<std::string, unsigned int> names, std::unordered_map<std::string, std::string> elements);

                /**
                 * @brief Get a value from the storage. 
                 *        Hydrogens will always return 0. 
                 *        If the key is not found, the average number of bonds in this residue for that element is returned. 
                 */
                double get(AtomKey key) const;

                /**
                 * @brief Insert a new element into the map. 
                 */
                void insert(AtomKey key, unsigned int value);

                /**
                 * @brief Insert a new element into the map. 
                 */
                void insert(std::string name, std::string symbol, unsigned int value);

            private: 
                std::unordered_map<AtomKey, unsigned int> map;      // the actual map data
                std::unordered_map<std::string, double> average;    // the average number of bonds for each element in this residue.
                
                void calculate_average(std::unordered_map<std::string, std::string> elements);
        };
    }
}