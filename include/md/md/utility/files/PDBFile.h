#pragma once

#include <io/File.h>

namespace md {
    // Protein data bank file
    struct PDBFile : public io::File {
        PDBFile() = default;
        PDBFile(const std::string& name) : File(name, "pdb") {}
        PDBFile(const char* name) : PDBFile(std::string(name)) {}
        ~PDBFile() override = default;
    };
}