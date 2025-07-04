// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <io/pdb/PDBStructure.h>

namespace ausaxs::io {
    struct Writer {
        static void write(const io::pdb::PDBStructure& s, const io::File&);
    };
}
