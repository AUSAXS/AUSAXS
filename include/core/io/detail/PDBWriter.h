#pragma once

#include <io/pdb/PDBFwd.h>
#include <io/IOFwd.h>

namespace ausaxs::io::detail::pdb {
    struct PDBWriter {
        static void write(const io::pdb::PDBStructure& data, const io::File& path);
    };
}