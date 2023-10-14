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
             * @brief Get the atomic group based on the residue name and atom type.
             */
            constants::atomic_group_t get_atomic_group(const std::string& name, constants::atom_t atom);

        private: 
            /**
             * @brief Insert a residue into the storage. 
             */
            void insert(const std::string& name, const detail::ResidueMap& residue);

            /**
             * @brief Initialize this storage. All residue files present in the storage directory will be loaded. 
             */
            void initialize();

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