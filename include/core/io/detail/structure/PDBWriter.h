// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <io/pdb/PDBFwd.h>

namespace ausaxs::io::detail::pdb {
    void write(const io::pdb::PDBStructure& data, const io::File& path);
}