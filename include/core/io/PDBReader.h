#pragma once

#include <io/pdb/PDBStructure.h>

namespace ausaxs::io::detail {
    class PDBReader {
        static io::pdb::PDBStructure read(const io::File& path);
    };
}