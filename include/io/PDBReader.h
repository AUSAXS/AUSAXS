#pragma once

class File;

#include "io/Reader.h"
#include "data/Terminate.h"
#include "data/Atom.h"
#include "data/Hetatom.h"

#include <fstream>

/**
 * @brief \class PDBReader. 
 * 
 * This class handles reading of input PDB format data files. 
 */
class PDBReader : public Reader {
public:
    /**
     * @brief Constructor.
     * @param file Path to the input PDB format data file. 
     */
    PDBReader(File* const file) : file(file) {}

    /**
     * @brief Read a PDB format data file.
     * @param input_path Path to the input PDB format data file. 
     */
    void read(const string& input_path) override;

private: 
    File* const file; // The File backing this Reader. 
};