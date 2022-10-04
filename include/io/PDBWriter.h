#pragma once

class File;

#include <io/Writer.h>
#include <data/Terminate.h>
#include <data/Atom.h>
#include <data/Water.h>

/**
 * @brief \class PDBWriter.
 *               This class handles writing a File object into a PDB format data file.
 */
class PDBWriter : public Writer {
    public:
        /**
         * @brief Constructor. 
         * @param file Path to where the backing File object will be saved. 
         */
        PDBWriter(File* file) : file(file) {}

        /**
         * @brief Write the backing File to disk. 
         * @param path the output path.
         */
        void write(std::string output_path) override;

    private: 
        File* file; // The File backing this Reader. 

        /**
         * @brief Create a string representation of this File.
         * @return The string representation. 
         */
        std::string as_pdb() const;
};