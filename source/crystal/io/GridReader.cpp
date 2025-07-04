// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <crystal/io/GridReader.h>
#include <utility/Exceptions.h>
#include <utility/Basis3D.h>

using namespace ausaxs;

std::pair<Basis3D, std::vector<Vector3<double>>> crystal::io::GridReader::read(const ::io::ExistingFile&) const {
    throw except::unexpected("GridReader::read: Not implemented");
}