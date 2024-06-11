/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/files/GROFile.h>
#include <md/utility/Exceptions.h>

#include <fstream>

using namespace md;

std::string GROFile::get_unit_cell() const {
    if (!exists()) {throw except::io_error("GROFile::get_unit_cell: \"" + path + "\" does not exist.");}
    std::ifstream in(path);

    // the unit cell is just the last line of the file
    std::string line, cell;
    while (std::getline(in, line)) {cell = line;}
    return cell;
}