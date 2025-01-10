/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <io/detail/PDBWriter.h>
#include <io/pdb/PDBStructure.h>
#include <io/File.h>
#include <io/pdb/Terminate.h>
#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <utility/Exceptions.h>
#include <settings/GeneralSettings.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::io::detail::pdb;
using namespace ausaxs::io::pdb;

std::vector<std::string> as_pdb(const PDBStructure& f) {
    std::vector<std::string> files;
    std::string s = f.header.get();
    Terminate t = f.terminate;

    unsigned int count = 0;
    int i_ter = t.serial;
    bool printed_ter = i_ter == -1;
    for (unsigned int i = 0; i < f.atoms.size(); i++) {
        if (static_cast<int>(i) == i_ter) { // check if this is where the terminate is supposed to go
            t.set_serial(t.serial % 100000);
            s += t.as_pdb(); // write it if so
            printed_ter = true;
        }
        s += f.atoms[i].as_pdb();
        count++;
        if (count == 100000) {
            count = 0;
            files.push_back(std::move(s));
            s = "";
        }
    }

    // print terminate if missing
    if (!printed_ter) {
        t.set_serial(t.serial % 100000);
        s += t.as_pdb();
    }

    // print hetatoms
    for (unsigned int i = 0; i < f.waters.size(); i++) {
        s += f.waters[i].as_pdb();
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

void PDBWriter::write(const PDBStructure& s, const io::File& path) {
    path.directory().create();

    auto content = as_pdb(s);
    if (content.size() == 1) {
        std::ofstream output(path);
        if (!output.is_open()) {throw std::ios_base::failure("PDBWriter::write: Could not open file \"" + path.str() + "\"");}
        output << content.at(0) << std::flush;
        output.close();
        if (settings::general::verbose) {std::cout << "Output written to file " + path.str() + "." << std::endl;}
    }
    else {
        for (unsigned int i = 0; i < content.size(); i++) {
            auto nfile = path.append("_part" + std::to_string(i+1));
            if (settings::general::verbose) {std::cout << "Output written to file " + nfile.str() << std::endl;}
            std::ofstream output(nfile);
            if (!output.is_open()) {throw std::ios_base::failure("PDBWriter::write: Could not open file \"" + path.str() + "\"");}
            output << content.at(i) << std::flush;
            output.close();
        }
    }
}