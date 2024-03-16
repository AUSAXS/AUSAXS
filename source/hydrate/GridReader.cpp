/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/GridReader.h>
#include <hydrate/GridMember.h>
#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/culling/CullingStrategy.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>
#include <io/ExistingFile.h>
#include <math/Vector3.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <utility/Limit3D.h>

#include <fstream>

using namespace hydrate;

grid::Grid GridReader::read(const io::ExistingFile& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw except::io_error("GridReader::read(const std::string& filename): Could not open file " + filename);
    }

    std::string line;
    std::getline(file, line);
    auto tokens = utility::split(line, ' ');
    Vector3 xaxis(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));

    std::getline(file, line);
    tokens = utility::split(line, ' ');
    Vector3 yaxis(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));

    std::getline(file, line);
    tokens = utility::split(line, ' ');
    Vector3 zaxis(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));

    // sanity checks
    if (xaxis.y() != 0 || xaxis.z() != 0) {
        throw except::io_error("GridReader::read(const std::string& filename): Grid x-axis must be parallel to the x-axis");
    }
    if (yaxis.x() != 0 || yaxis.z() != 0) {
        throw except::io_error("GridReader::read(const std::string& filename): Grid y-axis must be parallel to the y-axis");
    }
    if (zaxis.x() != 0 || zaxis.y() != 0) {
        throw except::io_error("GridReader::read(const std::string& filename): Grid z-axis must be parallel to the z-axis");
    }    
    if (xaxis.norm() == 0 || yaxis.norm() == 0 || zaxis.norm() == 0) {
        throw except::io_error("GridReader::read(const std::string& filename): Grid axes cannot be zero");
    }

    // Axis3D axes({0, xaxis.x(), 0, yaxis.y(), 0, zaxis.z()});
    // Grid grid(axes);
    // for (unsigned int x = 0; x < axes.x.max; ++x) {
    //     for (unsigned int y = 0; y < axes.y.max; ++y) {
    //         for (unsigned int z = 0; z < axes.z.max; ++z) {
    //             file >> grid.index(x, y, z);
    //         }
    //     }
    // }

    // return grid;
    return grid::Grid(Limit3D{1, 1, 1, 1, 1, 1});
}
