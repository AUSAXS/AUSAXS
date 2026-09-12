// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/Reader.h>

#include <io/detail/structure/CIFReader.h>
#include <io/detail/structure/PDBReader.h>
#include <io/detail/structure/XYZReader.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

io::pdb::PDBStructure io::Reader::read(const io::File& file) {
    auto ext = utility::to_lowercase(file.extension());
    if (file.extension() == ".xml" || file.extension() == ".XML") { // .xml PDBStructure
        throw except::invalid_argument("PDBStructure::construct_reader: .xml input PDBStructures are not supported.");
    }

    if (ext == ".pdb") { // .pdb structure
        return io::detail::pdb::read(file);
    }

    if (ext == ".ent") { // .ent structure
        return io::detail::pdb::read(file);
    }

    if (ext == ".cif") { // .cif structure
        return io::detail::cif::read(file);
    }

    if (ext == ".xyz") { // .xyz structure
        return io::detail::xyz::read(file);
    }

    // anything else - we cannot handle this
    if (auto format = constants::filetypes::detail::guess_type(file); !format.empty()) {
        throw except::invalid_argument("PDBStructure::construct_reader: Expected a structure file format, but got what appears to be a " + format + " file.");
    }
    throw except::invalid_argument("PDBStructure::construct_reader: Unsupported extension \"" + file.extension() + "\" of input file \"" + file.str() + "\".");
}