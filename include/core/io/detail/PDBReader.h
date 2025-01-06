#pragma once

#include <io/pdb/PDBStructure.h>

namespace ausaxs::io::detail::pdb {
    io::pdb::PDBStructure read(const io::File& path);
}