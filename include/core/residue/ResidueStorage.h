#pragma once

#include <residue/ResidueParser.h>

namespace residue {
    /**
     * @brief A storage container for residues.
     */
    class ResidueStorage {
        public: 
            /**
             * @brief Default constructor. The storage will be initialized with all residue files present in the storage directory.
             */
            ResidueStorage();

            ~ResidueStorage();

            /**
             * @brief Get a residue from the storage. If the residue is not found, it will be downloaded. 
             */
            detail::ResidueMap& get(const std::string& name);

            /**
             * @brief Check if a residue is present in the storage. 
             */
            bool contains(const std::string& name);

            /**
             * @brief Get the atomic group based on the residue name and atom type.
             * 
             * @param residue_name The name of the residue, e.g. GLY or ALA.
             * @param atom_name The name of the atom, e.g. CH2 or NH2
             * @param atom The atom type. This is required to avoid ambiguities since the name is always capitalized in PDB files, so otherwise we cannot distinguish between e.g. a C-alpha (CA) and a calcium (Ca).
             */
            constants::atomic_group_t get_atomic_group(const std::string& residue_name, const std::string& atom_name, constants::atom_t atom);

        private: 
            /**
             * @brief Insert a residue into the storage. 
             */
            void insert(const std::string& name, const detail::ResidueMap& residue);

            /**
             * @brief Initialize this storage. All residue files present in the storage directory will be loaded. 
             */
            void initialize();
            bool initialized = false;

            /**
             * @brief Download and load a residue from the web. 
             */
            void download_residue(const std::string& name);

            /**
             * @brief Write a residue to the auto-loaded file. 
             */
            void write_residue(const std::string& name);

            std::unordered_map<std::string, detail::ResidueMap> data;
    };
}