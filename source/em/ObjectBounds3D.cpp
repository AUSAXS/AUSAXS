// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <em/ObjectBounds3D.h>

#include <em/ObjectBounds2D.h>

#include <numeric>

using namespace ausaxs::em;

ObjectBounds3D::ObjectBounds3D(int size_x, int size_y, int size_z) : bounds(size_z, ObjectBounds2D(size_x, size_y)), _size_x(size_x), _size_y(size_y), _size_z(size_z) {}

ObjectBounds3D::~ObjectBounds3D() = default;

ObjectBounds2D& ObjectBounds3D::operator[](int z) {return bounds[z];}

const ObjectBounds2D& ObjectBounds3D::operator[](int z) const {return bounds[z];}

int ObjectBounds3D::total_volume() const {return _size_x*_size_y*_size_z;}

int ObjectBounds3D::bounded_volume() const {return std::accumulate(bounds.begin(), bounds.end(), 0, [] (int volume, const ObjectBounds2D& bound) {return volume += bound.bounded_area();});}

int ObjectBounds3D::size_x() const {return _size_x;}

int ObjectBounds3D::size_y() const {return _size_y;}

int ObjectBounds3D::size_z() const {return _size_z;}