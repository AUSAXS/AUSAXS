#include <em/ObjectBounds3D.h>
#include <em/ObjectBounds2D.h>

#include <numeric>

using namespace em;

ObjectBounds3D::ObjectBounds3D(unsigned int size_x, unsigned int size_y, unsigned int size_z) : bounds(size_z, ObjectBounds2D(size_x, size_y)), size_x(size_x), size_y(size_y), size_z(size_z) {}

ObjectBounds2D& ObjectBounds3D::operator[](unsigned int z) {return bounds[z];}

const ObjectBounds2D& ObjectBounds3D::operator[](unsigned int z) const {return bounds[z];}

unsigned int ObjectBounds3D::total_volume() const {return size_x*size_y*size_z;}

unsigned int ObjectBounds3D::bounded_volume() const {return std::accumulate(bounds.begin(), bounds.end(), 0, [] (unsigned int volume, const ObjectBounds2D& bound) {return volume += bound.bounded_area();});}
