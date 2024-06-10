#pragma once

#include <md/utility/files/File.h>

namespace gmx {
    // Protein data bank file
    struct PDBFile : public detail::File {
        PDBFile() = default;
        PDBFile(const std::string& name) : File(name, "pdb") {}
        PDBFile(const char* name) : PDBFile(std::string(name)) {}
    };
}