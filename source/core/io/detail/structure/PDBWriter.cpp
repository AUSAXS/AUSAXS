// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/detail/structure/PDBWriter.h>
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

    // at most this many records can be written to a single file before the serial wraps around and a new file must be started
    constexpr int serial_period = 100000;

    int n_serial = 0; // records written to the current file; all of them consume a serial
    auto next_file = [&] () {
        files.push_back(std::move(s));
        s = "";
        n_serial = 0;
    };

    int i_ter = t.serial;
    bool printed_ter = i_ter == -1;
    for (unsigned int i = 0; i < f.atoms.size(); i++) {
        if (static_cast<int>(i) == i_ter) { // check if this is where the terminate is supposed to go
            if (serial_period <= n_serial) {next_file();}
            t.set_serial(t.serial % serial_period);
            s += t.as_pdb(); // write it if so
            n_serial++;
            printed_ter = true;
        }

        if (serial_period <= n_serial) {next_file();}
        s += f.atoms[i].as_pdb();
        n_serial++;
    }

    // print terminate if missing
    if (!printed_ter) {
        if (serial_period <= n_serial) {next_file();}
        t.set_serial(t.serial % serial_period);
        s += t.as_pdb();
        n_serial++;
    }

    // print hetatoms
    for (unsigned int i = 0; i < f.waters.size(); i++) {
        if (serial_period <= n_serial) {next_file();}
        s += f.waters[i].as_pdb();
        n_serial++;
    }

    s += f.footer.get();
    files.push_back(std::move(s));
    return files;
}

void io::detail::pdb::write(const PDBStructure& s, const io::File& path) {
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