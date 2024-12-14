#pragma once

#include <io/File.h>
#include <io/pdb/PDBStructure.h>
#include <residue/ResidueFwd.h>

#include <vector>

namespace ausaxs::io::detail::cif {
    static io::pdb::PDBStructure read(const io::File& path);
    static std::vector<residue::detail::Residue> read_residue(const io::File& path);
}