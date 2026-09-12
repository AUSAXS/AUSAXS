// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/detail/structure/PDBWriter.h>

#include <io/File.h>
#include <io/pdb/PDBStructure.h>
#include <io/pdb/Terminate.h>
#include <settings/GeneralSettings.h>

#include <fstream>
#include <iostream>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::io::detail::pdb;
using namespace ausaxs::io::pdb;

namespace {
    std::vector<std::string> as_pdb(const PDBStructure& f) {
        std::vector<std::string> files;
        std::string s = f.header.get();
        Terminate t = f.terminate;

        int count = 0;
        int i_ter = t.serial;
        bool printed_ter = i_ter == -1;
        for (int i = 0; i < static_cast<int>(f.atoms.size()); i++) {
            if (i_ter == i) { // check if this is where the terminate is supposed to go
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
        for (const auto& water : f.waters) {
            s += water.as_pdb();
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
}

void io::detail::pdb::write(const PDBStructure& data, const io::File& path) {
    path.directory().create();

    auto content = as_pdb(data);
    if (content.size() == 1) {
        std::ofstream output(path);
        if (!output.is_open()) {throw ausaxs::except::io_error("PDBWriter::write: Could not open file \"" + path.str() + "\"");}
        output << content.at(0) << std::flush;
        output.close();
        if (settings::general::verbose) {std::cout << "Output written to file " + path.str() + "." << std::endl;}
    }
    else {
        for (int i = 0; i < static_cast<int>(content.size()); i++) {
            auto nfile = path.append("_part" + std::to_string(i+1));
            if (settings::general::verbose) {std::cout << "Output written to file " + nfile.str() << std::endl;}
            std::ofstream output(nfile);
            if (!output.is_open()) {throw ausaxs::except::io_error("PDBWriter::write: Could not open file \"" + path.str() + "\"");}
            output << content.at(i) << std::flush;
            output.close();
        }
    }
}