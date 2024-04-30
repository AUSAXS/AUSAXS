/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/ObjectBounds2D.h>
#include <utility/Limit.h>

#include <numeric>
#include <utility/Exceptions.h>

using namespace em;

ObjectBounds2D::ObjectBounds2D(unsigned int size_x, unsigned int size_y) : bounds(size_x, Limit(0, size_y)), N(size_x), M(size_y) {}

ObjectBounds2D::~ObjectBounds2D() = default;

void ObjectBounds2D::set_bounds(unsigned int x, const Limit& limit) {
    if (x >= N) {throw except::out_of_bounds("ObjectBounds2D::set_bounds: index " + std::to_string(x) + " is out of range (" + std::to_string(N) + ")");}
    if (limit.max > M) {throw except::out_of_bounds("ObjectBounds2D::set_bounds: limit " + std::to_string(limit.max) + " is out of bounds (" + std::to_string(M) + ")");}
    bounds[x] = limit;
}

void ObjectBounds2D::set_bounds(unsigned int x, unsigned int min, unsigned int max) {
    set_bounds(x, Limit(min, max));
}

void ObjectBounds2D::set_min(unsigned int x, unsigned int min) {
    if (x >= N) {throw except::out_of_bounds("ObjectBounds2D::set_min: index " + std::to_string(x) + " is out of range (" + std::to_string(N) + ")");}
    if (min > M) {throw except::out_of_bounds("ObjectBounds2D::set_min: limit " + std::to_string(min) + " is out of bounds (" + std::to_string(M) + ")");}
    bounds[x].min = min;
}

void ObjectBounds2D::set_max(unsigned int x, unsigned int max) {
    if (x >= N) {throw except::out_of_bounds("ObjectBounds2D::set_max: index " + std::to_string(x) + " is out of range (" + std::to_string(N) + ")");}
    if (max > M) {throw except::out_of_bounds("ObjectBounds2D::set_max: limit " + std::to_string(max) + " is out of bounds (" + std::to_string(M) + ")");}
    bounds[x].max = max;
}

const Limit& ObjectBounds2D::operator[](unsigned int x) const {return bounds[x];}

unsigned int ObjectBounds2D::size_x() const {return N;}

unsigned int ObjectBounds2D::size_y() const {return M;}

bool ObjectBounds2D::empty() const {return bounded_area() == N;}

unsigned int ObjectBounds2D::bounded_area() const {return std::accumulate(bounds.begin(), bounds.end(), 0u, [] (unsigned int area, const Limit& limit) {return area += static_cast<unsigned int>(limit.max+1 - limit.min);});}

unsigned int ObjectBounds2D::total_area() const {return N*M;}

bool ObjectBounds2D::operator==(const ObjectBounds2D& other) const = default;