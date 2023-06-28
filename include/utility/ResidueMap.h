#pragma once

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
            AtomKey(const std::string& name, const std::string& symbol);
            std::string name;
            std::string symbol;

            bool operator==(const AtomKey& other) const;
        };    
    }
}

namespace std {
    template <>
    struct hash<saxs::detail::AtomKey> {
        unsigned int operator()(const saxs::detail::AtomKey& k) const;
    };
}

namespace saxs {
    namespace detail {
        /**
         * @brief A simple map that stores the number of hydrogen bonds for each atom in a residue. 
         *        Since we do not model actual hydrogens, we instead modify the charge of the atom to which the hydrogen is bonded.
         *        This map serves as a lookup table for these effective charge calculations. 
         */
        class ResidueMap {
            public:
                /**
                 * @brief Default constructor.
                 */
                ResidueMap();

                /**
                 * @brief Construct a ResidueMap from an existing map.
                 */
                ResidueMap(const std::unordered_map<AtomKey, int>& map);

                /**
                 * @brief Get a value from the storage. 
                 *        Hydrogens will always return 0. 
                 *        If the key is not found, the average number of bonds in this residue for that element is returned. 
                 */
                double get(const AtomKey& key);

                /**
                 * @brief Get a value from the storage. 
                 *        Hydrogens will always return 0. 
                 *        If the key is not found, the average number of bonds in this residue for that element is returned. 
                 * 
                 * @param name The name of the atom.
                 * @param symbol The symbol of the atom.
                 */
                double get(const std::string& name, const std::string& symbol) {return this->get(AtomKey(name, symbol));}

                /**
                 * @brief Insert a new element into the map. 
                 * 
                 * @param key The key to insert.
                 * @param value The number of bonds.
                 */
                void insert(const AtomKey& key, int value);

                /**
                 * @brief Insert a new element into the map. 
                 * 
                 * @param name The name of the atom.
                 * @param symbol The symbol of the atom.
                 * @param value The number of bonds.
                 */
                void insert(const std::string& name, const std::string& symbol, int value);

                const std::unordered_map<AtomKey, int>& get_map() const;

                std::unordered_map<AtomKey, int>::const_iterator begin() const {return map.begin();}
                std::unordered_map<AtomKey, int>::const_iterator end() const {return map.end();}
                std::unordered_map<AtomKey, int>::iterator begin() {return map.begin();}
                std::unordered_map<AtomKey, int>::iterator end() {return map.end();}

            private: 
                std::unordered_map<AtomKey, int> map;               // the actual map data
                std::unordered_map<std::string, double> average;    // the average number of bonds for each element in this residue.
                bool update_average = false;                        // whether the average needs to be updated
                
                /**
                 * @brief Calculate the average number of bonds for this residue. 
                 */
                void calculate_average();
        };
    }
}