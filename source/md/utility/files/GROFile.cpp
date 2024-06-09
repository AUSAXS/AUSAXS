#include <utility/files/GROFile.h>
#include <utility/Exceptions.h>

#include <fstream>

using namespace gmx;

std::string GROFile::get_unit_cell() const {
    if (!exists()) {throw except::io_error("GROFile::get_unit_cell: \"" + path + "\" does not exist.");}
    std::ifstream in(path);

    // the unit cell is just the last line of the file
    std::string line, cell;
    while (std::getline(in, line)) {cell = line;}
    return cell;
}