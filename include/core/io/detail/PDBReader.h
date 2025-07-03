// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/pdb/PDBStructure.h>

namespace ausaxs::io::detail::pdb {
    io::pdb::PDBStructure read(const io::File& path);
}