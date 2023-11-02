#include <hist/distribution/Distribution1D.h>

#include <cmath>

using namespace hist;

Distribution1D::Distribution1D(unsigned int size) : data(size) {}
Distribution1D::Distribution1D(unsigned int size, type value) : data(size, value) {}

void Distribution1D::add(float distance, type value) {data.index(std::round(distance)) += value;}
void Distribution1D::add(int32_t i, type value) {data.index(i) += value;}

Distribution1D::type& Distribution1D::index(unsigned int i) {return data.index(i);}
const Distribution1D::type& Distribution1D::index(unsigned int i) const {return data.index(i);}

const typename std::vector<Distribution1D::type>::const_iterator Distribution1D::begin() const {return data.begin();}
const typename std::vector<Distribution1D::type>::const_iterator Distribution1D::end() const {return data.end();}

typename std::vector<Distribution1D::type>::iterator Distribution1D::begin() {return data.begin();}
typename std::vector<Distribution1D::type>::iterator Distribution1D::end() {return data.end();}

std::size_t Distribution1D::size() const {return data.size();}
bool Distribution1D::empty() const {return size() == 0;}

void Distribution1D::resize(unsigned int size) {data.resize(size);}

container::Container1D<Distribution1D::type>& Distribution1D::get_container() {return data;}