#include <hist/distribution/WeightedDistribution3D.h>

#include <cmath>

using namespace hist;

WeightedDistribution3D::WeightedDistribution3D(unsigned int size_x, unsigned int size_y, unsigned int size_z, type value) : data(size_x, size_y, size_z, value) {}

void WeightedDistribution3D::add(unsigned int x, unsigned int y, float distance, type value) {
    data.index(x, y, std::round(distance)).add(distance, value);
}

WeightedDistribution3D::type& WeightedDistribution3D::index(unsigned int x, unsigned int y, unsigned int z) {return data.index(x, y, z).value;}
const WeightedDistribution3D::type& WeightedDistribution3D::index(unsigned int x, unsigned int y, unsigned int z) const {return data.index(x, y, z).value;}

const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution3D::begin(unsigned int x, unsigned int y) const {return data.begin(x, y);}
const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution3D::end(unsigned int x, unsigned int y) const {return data.end(x, y);}
const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution3D::begin() const {return data.begin();}
const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution3D::end() const {return data.end();}

typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution3D::begin(unsigned int x, unsigned int y) {return data.begin(x, y);}
typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution3D::end(unsigned int x, unsigned int y) {return data.end(x, y);}
typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution3D::begin() {return data.begin();}
typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution3D::end() {return data.end();}

std::size_t WeightedDistribution3D::size_x() const {return data.size_x();}
std::size_t WeightedDistribution3D::size_y() const {return data.size_y();}
std::size_t WeightedDistribution3D::size_z() const {return data.size_z();}

void WeightedDistribution3D::resize(unsigned int size) {data.resize(size);}