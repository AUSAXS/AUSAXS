#pragma once

#include <io/pdb/PDBStructure.h>

namespace ausaxs::io::detail::pdb {
    static io::pdb::PDBStructure read(const io::File& path);
}