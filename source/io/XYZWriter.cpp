/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/XYZWriter.h>
#include <data/Molecule.h>
#include <io/File.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <data/record/Atom.h>

#include <iomanip>

io::XYZWriter::XYZWriter(const io::File& path) : path(path) {
    path.directory().create();
    file.open(path);
    if (!file.is_open()) {
        throw except::io_error("io::XYZWriter: Could not open file " + path);
    }
}

io::XYZWriter::~XYZWriter() {
    file.close();
    console::print_info("Trajectory written to " + path);
}

void io::XYZWriter::write_frame(const data::Molecule* protein) {
    static unsigned int frame = 0;
    auto atoms = protein->get_atoms();
    file << " " << atoms.size() << std::endl;
    file << " Frame " << frame++ << std::endl;
    for (const auto& atom : atoms) {
        file << std::setw(10) << atom.serial << " " 
                << std::setw(10) << std::setprecision(6) << atom.coords.x() << " " 
                << std::setw(10) << std::setprecision(6) << atom.coords.y() << " " 
                << std::setw(10) << std::setprecision(6) << atom.coords.z() << std::endl;
    }
}