#pragma once

#include <io/pdb/PDBFwd.h>
#include <io/IOFwd.h>

namespace ausaxs::io::detail {
    struct PDBWriter {
        static void write(const io::pdb::PDBStructure&, const io::File& path);
    };
}