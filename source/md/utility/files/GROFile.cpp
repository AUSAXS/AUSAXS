/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/files/GROFile.h>
#include <md/utility/Exceptions.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::md;

std::string GROFile::get_unit_cell() const {
    if (!exists()) {throw except::io_error("GROFile::get_unit_cell: \"" + path() + "\" does not exist.");}
    std::ifstream in(path());

    // the unit cell is just the last line of the file
    std::string line, cell;
    while (std::getline(in, line)) {cell = line;}
    return cell;
}

unsigned int GROFile::size_solvent() const {
    if (!exists()) {throw except::io_error("GROFile::size_solvent: \"" + path() + "\" does not exist.");}
    std::ifstream in(path());

    std::string line;
    std::getline(in, line);
    std::getline(in, line);
    int total = std::stoi(line);

    // skip the coordinates
    int count_mol = 0;
    while(std::getline(in, line)) {
        if (line.substr(5, 3) == "SOL") {break;}
        ++count_mol;
    }

    return (total - count_mol)/4;
}