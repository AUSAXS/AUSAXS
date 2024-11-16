#pragma once

#include <constants/ConstantsFwd.h>

#include <unordered_map>
#include <string>

namespace ausaxs::residue::detail {
    /**
     * @brief The key type used in the SimpleResidueMap. 
     *        The atom type is required to avoid ambiguities since the name is always capitalized in PDB files,
     *        so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca). 
     */
    struct AtomKey {
        AtomKey(const std::string& name, constants::atom_t atom);
        std::string name;
        constants::atom_t atom;

        bool operator==(const AtomKey& other) const;
    };    
}

namespace std {
    template <>
    struct hash<ausaxs::residue::detail::AtomKey> {
        unsigned int operator()(const ausaxs::residue::detail::AtomKey& k) const;
    };
}

namespace ausaxs::residue::detail {
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
             * @brief Get the number of bound hydrogens.
             *        Hydrogens will always return 0. 
             *        If the key is not found, the average number of bonds in this residue for that element is returned. 
             * 
             * @param atom_name The name of the atom, e.g. CH2 or NH2
             * @param symbol The atom type. This is required to avoid ambiguities since the name is always capitalized in PDB files, so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca).
             */
            double get(const std::string& atom_name, constants::atom_t atom) {return this->get(AtomKey(atom_name, atom));}

            /**
             * @brief Check if a key is present in the map. 
             */
            bool contains(const std::string& name, constants::atom_t atom) const;

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
             * @param atom_name The name of the atom, e.g. CH2 or NH2
             * @param atom The atom type. This is required to avoid ambiguities since the name is always capitalized in PDB files, so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca).
             * @param value The number of bonds.
             */
            void insert(const std::string& atom_name, constants::atom_t atom, int value);

            const std::unordered_map<AtomKey, int>& get_map() const;

            /**
             * @brief Get the atomic group a given atom belongs to. 
             * 
             * @param atom_name The atom name, e.g. CH2 or NH2
             * @param atom_type The atom type. This is required to avoid ambiguities since the name is always capitalized in PDB files, so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca).
             */
            constants::atomic_group_t get_atomic_group(const std::string& atom_name, constants::atom_t atom_type);

            std::unordered_map<AtomKey, int>::const_iterator begin() const {return map.begin();}
            std::unordered_map<AtomKey, int>::const_iterator end() const {return map.end();}
            std::unordered_map<AtomKey, int>::iterator begin() {return map.begin();}
            std::unordered_map<AtomKey, int>::iterator end() {return map.end();}

        private: 
            std::unordered_map<AtomKey, int> map;                   // the actual map data
            std::unordered_map<constants::atom_t, double> average;  // the average number of bonds for each element in this residue.
            bool update_average = false;                            // whether the average needs to be updated
            
            /**
             * @brief Calculate the average number of bonds for this residue. 
             */
            void calculate_average();
    };
}