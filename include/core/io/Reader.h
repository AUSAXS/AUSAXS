#pragma once

#include <io/IOFwd.h>
#include <io/pdb/PDBStructure.h>

namespace ausaxs::io::detail {
    struct Reader {
        static io::pdb::PDBStructure read(const io::ExistingFile& file);
    };
}