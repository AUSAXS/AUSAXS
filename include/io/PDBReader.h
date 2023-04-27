#pragma once

#include <io/Reader.h>
#include <data/Terminate.h>
#include <data/Atom.h>
#include <data/Water.h>

class ProteinFile;

/**
 * @brief This class handles reading of input PDB format data files. 
 */
class PDBReader : public Reader {
    public:
        /**
         * @brief Constructor.
         * @param file Path to the input PDB format data file. 
         */
        PDBReader(ProteinFile* const file) : file(file) {}

        ~PDBReader() override = default;

        /**
         * @brief Read a PDB format data file.
         * 
         * @param path Path to the input PDB format data file. 
         */
        void read(const io::ExistingFile& path) override;

    private: 
        ProteinFile* const file; // The File backing this Reader. 
};