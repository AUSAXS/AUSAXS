// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/File.h>
#include <io/pdb/PDBStructure.h>
#include <residue/ResidueFwd.h>

#include <vector>

namespace ausaxs::io::detail::cif {
    io::pdb::PDBStructure read(const io::File& path);
    std::vector<residue::detail::Residue> read_residue(const io::File& path);
}