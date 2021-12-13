#pragma once

class File;

#include "io/Reader.h"
#include "data/Terminate.h"
#include "data/Atom.h"
#include "data/Hetatom.h"

#include <fstream>

class PDBReader : public Reader {
public:
    PDBReader(File* const file) : file(file) {}

    void read(const string& input_path) override;

private: 
    File* const file;
};