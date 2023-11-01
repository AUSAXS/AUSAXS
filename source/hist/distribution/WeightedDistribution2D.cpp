#include <hist/distribution/WeightedDistribution2D.h>

#include <cmath>

using namespace hist;

WeightedDistribution2D::WeightedDistribution2D(unsigned int size_x, unsigned int size_y, type value) : data(size_x, size_y, value) {}

void WeightedDistribution2D::add(unsigned int x, float distance, type value) {
    data.index(x, std::round(distance)).add(distance, value);
}

WeightedDistribution2D::type& WeightedDistribution2D::index(unsigned int x, unsigned int y) {return data.index(x, y).value;}
const WeightedDistribution2D::type& WeightedDistribution2D::index(unsigned int x, unsigned int y) const {return data.index(x, y).value;}

const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution2D::begin(unsigned int x) const {return data.begin(x);}
const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution2D::end(unsigned int x) const {return data.end(x);}
const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution2D::begin() const {return data.begin();}
const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution2D::end() const {return data.end();}

typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution2D::begin(unsigned int x) {return data.begin(x);}
typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution2D::end(unsigned int x) {return data.end(x);}
typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution2D::begin() {return data.begin();}
typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution2D::end() {return data.end();}

std::size_t WeightedDistribution2D::size_x() const {return data.size_x();}
std::size_t WeightedDistribution2D::size_y() const {return data.size_y();}

void WeightedDistribution2D::resize(unsigned int size) {data.resize(size);}