#pragma once

#include <io/IOFwd.h>
#include <data/Molecule.h>

namespace ausaxs::io::detail {
    class Reader {
        static data::Molecule read(const io::ExistingFile& file);
    };
}