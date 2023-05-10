#include <crystal/io/GridReader.h>
#include <utility/Exceptions.h>
#include <utility/Basis3D.h>

std::pair<Basis3D, std::vector<Vector3<double>>> crystal::io::GridReader::read(const ::io::ExistingFile&) const {
    throw except::unexpected("GridReader::read: Not implemented");
}