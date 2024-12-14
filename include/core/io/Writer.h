#pragma once

#include <io/IOFwd.h>
#include <io/pdb/PDBStructure.h>

namespace ausaxs::io::detail {
    struct Writer {
        virtual void write(const io::pdb::PDBStructure& s, const io::File&);
    };
}
