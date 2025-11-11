// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/Reader.h>
#include <io/detail/structure/CIFReader.h>
#include <io/detail/structure/PDBReader.h>
#include <io/detail/structure/XYZReader.h>
#include <io/ExistingFile.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

io::pdb::PDBStructure io::Reader::read(const io::File& path) {
    auto ext = utility::to_lowercase(path.extension());
    if (path.extension() == ".xml" || path.extension() == ".XML") { // .xml PDBStructure
        throw except::invalid_argument("PDBStructure::construct_reader: .xml input PDBStructures are not supported.");
    } else if (ext == ".pdb") { // .pdb structure
        return io::detail::pdb::read(path);
    } else if (ext == ".ent") { // .ent structure
        return io::detail::pdb::read(path);
    } else if (ext == ".cif") { // .cif structure
        return io::detail::cif::read(path);
    } else if (ext == ".xyz") { // .xyz structure
        return io::detail::xyz::read(path);
    } else { // anything else - we cannot handle this
        throw except::invalid_argument("PDBStructure::construct_reader: Unsupported extension \"" + path.extension() + "\" of input file \"" + path.str() + "\".");
    }
}