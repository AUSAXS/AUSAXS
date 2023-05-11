#include <io/XYZWriter.h>
#include <data/Protein.h>
#include <io/File.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <data/Atom.h>

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

void io::XYZWriter::write_frame(const Protein* protein) {
    static unsigned int frame = 0;
    auto atoms = protein->atoms();
    file << " " << atoms.size() << std::endl;
    file << " Frame " << frame++ << std::endl;
    for (const auto& atom : atoms) {
        file << std::setw(10) << atom.serial << " " 
                << std::setw(10) << std::setprecision(6) << atom.coords.x() << " " 
                << std::setw(10) << std::setprecision(6) << atom.coords.y() << " " 
                << std::setw(10) << std::setprecision(6) << atom.coords.z() << std::endl;
    }
}