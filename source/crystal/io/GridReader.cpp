#include <crystal/io/GridReader.h>

#include <fstream>

std::pair<Basis3D, std::vector<Vector3<double>>> crystal::io::GridReader::read(const std::string&) const {
    throw except::unexpected("GridReader::read: Not implemented");
}