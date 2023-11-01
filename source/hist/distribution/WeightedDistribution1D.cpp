#include <hist/distribution/WeightedDistribution1D.h>

#include <cmath>

using namespace hist;

WeightedDistribution1D::WeightedDistribution1D(unsigned int size, type value) : data(size, value) {}

void WeightedDistribution1D::add(float distance, type value) {data.index(std::round(distance)).add(distance, value);}

WeightedDistribution1D::type& WeightedDistribution1D::index(unsigned int i) {return data.index(i).value;}
const WeightedDistribution1D::type& WeightedDistribution1D::index(unsigned int i) const {return data.index(i).value;}

const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution1D::begin() const {return data.begin();}
const typename std::vector<hist::detail::WeightedEntry>::const_iterator WeightedDistribution1D::end() const {return data.end();}

typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution1D::begin() {return data.begin();}
typename std::vector<hist::detail::WeightedEntry>::iterator WeightedDistribution1D::end() {return data.end();}

std::size_t WeightedDistribution1D::size() const {return data.size();}
void WeightedDistribution1D::resize(unsigned int size) {data.resize(size);}

std::vector<WeightedDistribution1D::type>& WeightedDistribution1D::get_counts() {
    static std::vector<type> result;
    result = std::vector<type>(size());
    for (unsigned int i = 0; i < size(); ++i) {
        result[i] = data.index(i).value;
    }
    return result;
}