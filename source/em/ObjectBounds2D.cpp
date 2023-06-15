#include <em/ObjectBounds2D.h>
#include <utility/Limit.h>

#include <numeric>

using namespace em;

ObjectBounds2D::ObjectBounds2D(unsigned int size_x, unsigned int size_y) : bounds(size_x, Limit(0, size_y)), N(size_x), M(size_y) {}

ObjectBounds2D::~ObjectBounds2D() = default;

Limit& ObjectBounds2D::operator[](unsigned int x) {return bounds[x];}

const Limit& ObjectBounds2D::operator[](unsigned int x) const {return bounds[x];}

unsigned int ObjectBounds2D::size() const {return bounds.size();}

bool ObjectBounds2D::empty() const {return bounds.empty();}

unsigned int ObjectBounds2D::bounded_area() const {return std::accumulate(bounds.begin(), bounds.end(), 0, [] (unsigned int area, const Limit& limit) {return area += limit.max+1 - limit.min;});}

unsigned int ObjectBounds2D::total_area() const {return N*M;}

bool ObjectBounds2D::operator==(const ObjectBounds2D& other) const = default;