// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <io/ExistingFile.h>
#include <io/pdb/PDBStructure.h>

namespace ausaxs::io {
    struct Reader {
        static io::pdb::PDBStructure read(const io::File& file);
    };
}