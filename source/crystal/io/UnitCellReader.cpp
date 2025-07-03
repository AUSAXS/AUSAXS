// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <crystal/io/UnitCellReader.h>
#include <utility/StringUtils.h>
#include <io/ExistingFile.h>
#include <math/Vector3.h>
#include <utility/Basis3D.h>
#include <constants/Constants.h>

#include <fstream>

using namespace ausaxs;

std::pair<Basis3D, std::vector<Vector3<double>>> crystal::io::UnitCellReader::read(const ausaxs::io::ExistingFile& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) {throw except::io_error("GridReader::read: Could not open file " + filename.str());}

    std::string line;
    std::getline(file, line);
    if (line.substr(0, 6) != "BASIS") {throw except::io_error("GridReader::read: File \"" + filename.str() + "\" is not a crystal file (missing section: \"BASIS\").");}

    std::getline(file, line);
    auto tokens = utility::split(line, ' ');
    Vector3 xaxis(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));

    std::getline(file, line);
    tokens = utility::split(line, ' ');
    Vector3 yaxis(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));

    std::getline(file, line);
    tokens = utility::split(line, ' ');
    Vector3 zaxis(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));

    std::getline(file, line);
    std::getline(file, line);
    if (line.substr(0, 11) != "CRYSTALDATA") {throw except::io_error("GridReader::read: File \"" + filename.str() + "\" is not a crystal file (missing section: \"CRYSTALDATA\").");}

    // sanity checks
    if (xaxis.y() != 0 || xaxis.z() != 0) {throw except::io_error("GridReader::read: Grid x-axis must be parallel to the x-axis");}
    if (yaxis.x() != 0 || yaxis.z() != 0) {throw except::io_error("GridReader::read: Grid y-axis must be parallel to the y-axis");}
    if (zaxis.x() != 0 || zaxis.y() != 0) {throw except::io_error("GridReader::read: Grid z-axis must be parallel to the z-axis");}
    if (xaxis.norm() == 0 || yaxis.norm() == 0 || zaxis.norm() == 0) {throw except::io_error("GridReader::read: Grid axes cannot be zero");}

    // calculate the distance between the edges of the box
    double distance = std::sqrt(std::pow(xaxis.x(), 2) + std::pow(yaxis.y(), 2) + std::pow(zaxis.z(), 2));
    if (distance > constants::axes::d_axis.max) {
        throw except::io_error("PDBReader::read: The distance between the edges of the box is " + std::to_string(distance) + " Å, which is larger than the maximum allowed distance of " + std::to_string(constants::axes::d_axis.max) + " Å.");
    }

    std::getline(file, line); // skip empty line
    std::vector<Vector3<double>> voxels;
    while (std::getline(file, line)) {
        tokens = utility::split(line, ' ');
        voxels.push_back(Vector3<double>(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2])));
    }

    return std::make_pair(Basis3D(xaxis, yaxis, zaxis), std::move(voxels));
}