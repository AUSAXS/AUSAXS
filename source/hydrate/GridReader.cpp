#include <hydrate/GridReader.h>
#include <utility/Utility.h>

#include <fstream>

using namespace hydrate;

Grid GridReader::read(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("GridReader::read(const std::string& filename): Could not open file " + filename);
    }

    std::string line;
    std::getline(file, line);
    auto tokens = utility::split(line, ' ');
    double xdim = std::stod(tokens[0]);
    double ydim = std::stod(tokens[1]);
    double zdim = std::stod(tokens[2]);

    // Grid grid(xdim, ydim, zdim);
    // for (unsigned int x = 0; x < xdim; ++x) {
    //     for (unsigned int y = 0; y < ydim; ++y) {
    //         for (unsigned int z = 0; z < zdim; ++z) {
    //             file >> grid.index(x, y, z);
    //         }
    //     }
    // }

    // return grid;
    return Grid({1, 1, 1, 1, 1, 1}, 1);
}
