#pragma once

class File;

#include "io/Writer.h"
#include "data/Terminate.h"
#include "data/Atom.h"
#include "data/Hetatom.h"

#include <fstream>

class PDBWriter : public Writer {
public:
    PDBWriter(File* const file) : file(file) {}

    /**
     * @brief write this File to disk. 
     * @param path the output path.
     */
    void write(const string& output_path) override;

private: 
    File* const file;

    /**
     * @brief Create a string representation of this File.
     * @return The string representation. 
     */
    string as_pdb() const;
};