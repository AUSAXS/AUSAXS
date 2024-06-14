/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/PDBWriter.h>
#include <io/Writer.h>
#include <data/detail/AtomCollection.h>
#include <io/File.h>
#include <data/record/Terminate.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <utility/Exceptions.h>
#include <settings/GeneralSettings.h>

#include <fstream>

using namespace io::detail;
using namespace data::record;

PDBWriter::PDBWriter(data::detail::AtomCollection* file) : file(file) {}

PDBWriter::~PDBWriter() = default;

void PDBWriter::write(const io::File& path) {
    path.directory().create();

    auto content = as_pdb();
    if (content.size() == 1) {
        std::ofstream output(path);
        if (!output.is_open()) {throw std::ios_base::failure("PDBWriter::write: Could not open file \"" + path + "\"");}
        output << content[0] << std::flush;
        output.close();
        if (settings::general::verbose) {std::cout << "Output written to file " + path + "." << std::endl;}
    }
    else {
        for (unsigned int i = 0; i < content.size(); i++) {
            auto nfile = path.append("_part" + std::to_string(i+1));
            if (settings::general::verbose) {std::cout << "Output written to file " + nfile << std::endl;}
            std::ofstream output(nfile);
            if (!output.is_open()) {throw std::ios_base::failure("PDBWriter::write: Could not open file \"" + path + "\"");}
            output << content[i] << std::flush;
            output.close();
        }
    }
}

std::vector<std::string> PDBWriter::as_pdb() const {
    data::detail::AtomCollection& f = *file;
    std::vector<std::string> files;
    std::string s = f.header.get();

    unsigned int count = 0;
    unsigned int i_ter = f.terminate.serial;
    bool printed_ter = false;
    for (unsigned int i = 0; i < f.protein_atoms.size(); i++) {
        if (i == i_ter) { // check if this is where the terminate is supposed to go
            f.terminate.set_serial(f.terminate.serial % 100000);
            s += f.terminate.as_pdb(); // write it if so
            printed_ter = true;
        }
        s += f.protein_atoms[i].as_pdb();
        count++;
        if (count == 100000) {
            count = 0;
            files.push_back(std::move(s));
            s = "";
        }
    }

    // print terminate if missing
    if (!printed_ter) {
        f.terminate.set_serial(f.terminate.serial % 100000);
        s += f.terminate.as_pdb();
    }

    // print hetatoms
    for (unsigned int i = 0; i < f.hydration_atoms.size(); i++) {
        s += f.hydration_atoms[i].as_pdb();
        count++;
        if (count == 100000) {
            count = 0;
            files.push_back(std::move(s));
            s = "";
        }
    }

    s += f.footer.get();
    files.push_back(std::move(s));
    return files;
}