#include <crystal/io/UnitCellReader.h>
#include <utility/Utility.h>

#include <fstream>

std::pair<Basis3D, std::vector<Vector3<double>>> crystal::io::UnitCellReader::read(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) {throw except::io_error("GridReader::read: Could not open file " + filename);}

    std::string line;
    std::getline(file, line);
    if (line.substr(0, 6) != "BASIS") {throw except::io_error("GridReader::read: File \"" + filename + "\" is not a crystal file (missing section: \"BASIS\").");}

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
    if (line.substr(0, 11) != "CRYSTALDATA") {throw except::io_error("GridReader::read: File \"" + filename + "\" is not a crystal file (missing section: \"CRYSTALDATA\").");}

    // sanity checks
    if (xaxis.y() != 0 || xaxis.z() != 0) {throw except::io_error("GridReader::read: Grid x-axis must be parallel to the x-axis");}
    if (yaxis.x() != 0 || yaxis.z() != 0) {throw except::io_error("GridReader::read: Grid y-axis must be parallel to the y-axis");}
    if (zaxis.x() != 0 || zaxis.y() != 0) {throw except::io_error("GridReader::read: Grid z-axis must be parallel to the z-axis");}
    if (xaxis.norm() == 0 || yaxis.norm() == 0 || zaxis.norm() == 0) {throw except::io_error("GridReader::read: Grid axes cannot be zero");}

    std::getline(file, line); // skip empty line
    std::vector<Vector3<double>> voxels;
    while (std::getline(file, line)) {
        tokens = utility::split(line, ' ');
        voxels.push_back(Vector3<double>(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2])));
    }

    return std::make_pair(Basis3D(xaxis, yaxis, zaxis), std::move(voxels));
}