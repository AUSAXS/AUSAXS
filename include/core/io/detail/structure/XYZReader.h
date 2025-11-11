// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <io/pdb/PDBStructure.h>

namespace ausaxs::io::detail::xyz {
    io::pdb::PDBStructure read(const io::File& path);
}