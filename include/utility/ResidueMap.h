#pragma once

#include <utility/Utility.h>

#include <unordered_map>
#include <string>

namespace saxs {
    namespace detail {
        /**
         * @brief The key type used in the SimpleResidueMap. 
         *        The symbol is required to avoid ambiguities since the name is always capitalized in PDB files,
         *        so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca). 
         */
        struct AtomKey {
            AtomKey(std::string name, std::string symbol) : name(utility::to_lowercase(name)), symbol(symbol) {}
            std::string name;
            std::string symbol;

            bool operator==(const AtomKey& other) const {
                return name == other.name;
            }
        };    
    }
}

namespace std {
    template <>
    struct hash<saxs::detail::AtomKey> {
        std::size_t operator()(const saxs::detail::AtomKey& k) const {return std::hash<std::string>()(k.name);}
    };
}

namespace saxs {
    namespace detail {
        class ResidueMap {
            public:
                /**
                 * @brief Create a new empty SimpleResidueMap.
                 */
                ResidueMap();

                ResidueMap(std::unordered_map<AtomKey, unsigned int> map);

                /**
                 * @brief Get a value from the storage. 
                 *        Hydrogens will always return 0. 
                 *        If the key is not found, the average number of bonds in this residue for that element is returned. 
                 */
                double get(AtomKey key);

                /**
                 * @brief Get a value from the storage. 
                 *        Hydrogens will always return 0. 
                 *        If the key is not found, the average number of bonds in this residue for that element is returned. 
                 */
                double get(std::string name, std::string symbol) {return this->get(AtomKey(name, symbol));}

                /**
                 * @brief Insert a new element into the map. 
                 */
                void insert(AtomKey key, unsigned int value);

                /**
                 * @brief Insert a new element into the map. 
                 */
                void insert(std::string name, std::string symbol, unsigned int value);

                std::unordered_map<AtomKey, unsigned int>::const_iterator begin() const {return map.begin();}
                std::unordered_map<AtomKey, unsigned int>::const_iterator end() const {return map.end();}
                std::unordered_map<AtomKey, unsigned int>::iterator begin() {return map.begin();}
                std::unordered_map<AtomKey, unsigned int>::iterator end() {return map.end();}

            private: 
                std::unordered_map<AtomKey, unsigned int> map;      // the actual map data
                std::unordered_map<std::string, double> average;    // the average number of bonds for each element in this residue.
                bool update_average = false;                        // whether the average needs to be updated
                
                void calculate_average();
        };
    }
}