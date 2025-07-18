// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/detail/XYZWriter.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/File.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>

#include <iomanip>

using namespace ausaxs::io::detail::xyz;

XYZWriter::XYZWriter(const io::File& path) : path(path) {
    path.directory().create();
    file.open(path);
    if (!file.is_open()) {
        throw except::io_error("io::XYZWriter: Could not open file " + path.str());
    }
}

XYZWriter::~XYZWriter() {
    file.close();
    console::print_info("Trajectory written to " + path);
}

void XYZWriter::write_frame(const data::Molecule* protein) {
    static unsigned int frame = 0;
    std::vector<data::AtomFF> atoms;
    atoms.reserve(protein->size_atom());
    for (const auto& body : protein->get_bodies()) {
        auto bsym = body.symmetry().explicit_structure();
        atoms.insert(atoms.end(), bsym.atoms.begin(), bsym.atoms.end());
    }
    file << " " << atoms.size() << std::endl;
    file << " Frame " << frame++ << std::endl;
    int i = 0;
    for (const auto& atom : atoms) {
        file << std::setw(10) << i++ << " " 
                << std::setw(10) << std::setprecision(6) << atom.coordinates().x() << " " 
                << std::setw(10) << std::setprecision(6) << atom.coordinates().y() << " " 
                << std::setw(10) << std::setprecision(6) << atom.coordinates().z() << std::endl;
    }
}