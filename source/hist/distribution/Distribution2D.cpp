#include <hist/distribution/Distribution2D.h>

#include <cmath>

using namespace hist;

Distribution2D::Distribution2D(unsigned int size_x, unsigned int size_y, type value) : data(size_x, size_y, value) {}

void Distribution2D::add(unsigned int x, float distance, type value) {data.index(x, std::round(distance)) += value;}
void Distribution2D::add(unsigned int x, int32_t i, type value) {data.index(x, i) += value;}

Distribution2D::type& Distribution2D::index(unsigned int x, unsigned int y) {return data.index(x, y);}
const Distribution2D::type& Distribution2D::index(unsigned int x, unsigned int y) const {return data.index(x, y);}

const typename std::vector<Distribution2D::type>::const_iterator Distribution2D::begin(unsigned int x) const {return data.begin(x);}
const typename std::vector<Distribution2D::type>::const_iterator Distribution2D::end(unsigned int x) const {return data.end(x);}
const typename std::vector<Distribution2D::type>::const_iterator Distribution2D::begin() const {return data.begin();}
const typename std::vector<Distribution2D::type>::const_iterator Distribution2D::end() const {return data.end();}

typename std::vector<Distribution2D::type>::iterator Distribution2D::begin(unsigned int x) {return data.begin(x);}
typename std::vector<Distribution2D::type>::iterator Distribution2D::end(unsigned int x) {return data.end(x);}
typename std::vector<Distribution2D::type>::iterator Distribution2D::begin() {return data.begin();}
typename std::vector<Distribution2D::type>::iterator Distribution2D::end() {return data.end();}

std::size_t Distribution2D::size_x() const {return data.size_x();}
std::size_t Distribution2D::size_y() const {return data.size_y();}
bool Distribution2D::empty() const {return size_x() == 0 || size_y() == 0;}

void Distribution2D::resize(unsigned int size) {data.resize(size);}