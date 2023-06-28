#include <em/ObjectBounds3D.h>
#include <em/ObjectBounds2D.h>
#include <utility/Limit.h>

#include <numeric>

using namespace em;

ObjectBounds3D::ObjectBounds3D(unsigned int size_x, unsigned int size_y, unsigned int size_z) : bounds(size_z, ObjectBounds2D(size_x, size_y)), _size_x(size_x), _size_y(size_y), _size_z(size_z) {}

ObjectBounds3D::~ObjectBounds3D() = default;

ObjectBounds2D& ObjectBounds3D::operator[](unsigned int z) {return bounds[z];}

const ObjectBounds2D& ObjectBounds3D::operator[](unsigned int z) const {return bounds[z];}

unsigned int ObjectBounds3D::total_volume() const {return _size_x*_size_y*_size_z;}

unsigned int ObjectBounds3D::bounded_volume() const {return std::accumulate(bounds.begin(), bounds.end(), 0, [] (unsigned int volume, const ObjectBounds2D& bound) {return volume += bound.bounded_area();});}

unsigned int ObjectBounds3D::size_x() const {return _size_x;}

unsigned int ObjectBounds3D::size_y() const {return _size_y;}

unsigned int ObjectBounds3D::size_z() const {return _size_z;}