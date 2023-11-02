#include <hist/distribution/Distribution3D.h>

#include <cmath>

using namespace hist;

Distribution3D::Distribution3D(unsigned int size_x, unsigned int size_y, unsigned int size_z, Distribution3D::type value) : data(size_x, size_y, size_z, value) {}

void Distribution3D::add(unsigned int x, unsigned int y, float distance, Distribution3D::type value) {data.index(x, y, std::round(distance)) += value;}
void Distribution3D::add(unsigned int x, unsigned int y, int32_t i, Distribution3D::type value) {data.index(x, y, i) += value;}

Distribution3D::type& Distribution3D::index(unsigned int x, unsigned int y, unsigned int z) {return data.index(x, y, z);}
const Distribution3D::type& Distribution3D::index(unsigned int x, unsigned int y, unsigned int z) const {return data.index(x, y, z);}

const typename std::vector<Distribution3D::type>::const_iterator Distribution3D::begin(unsigned int x, unsigned int y) const {return data.begin(x, y);}
const typename std::vector<Distribution3D::type>::const_iterator Distribution3D::end(unsigned int x, unsigned int y) const {return data.end(x, y);}
const typename std::vector<Distribution3D::type>::const_iterator Distribution3D::begin() const {return data.begin();}
const typename std::vector<Distribution3D::type>::const_iterator Distribution3D::end() const {return data.end();}

typename std::vector<Distribution3D::type>::iterator Distribution3D::begin(unsigned int x, unsigned int y) {return data.begin(x, y);}
typename std::vector<Distribution3D::type>::iterator Distribution3D::end(unsigned int x, unsigned int y) {return data.end(x, y);}
typename std::vector<Distribution3D::type>::iterator Distribution3D::begin() {return data.begin();}
typename std::vector<Distribution3D::type>::iterator Distribution3D::end() {return data.end();}

std::size_t Distribution3D::size_x() const {return data.size_x();}
std::size_t Distribution3D::size_y() const {return data.size_y();}
std::size_t Distribution3D::size_z() const {return data.size_z();}
bool Distribution3D::empty() const {return size_x() == 0 || size_y() == 0 || size_z() == 0;}

void Distribution3D::resize(unsigned int size) {data.resize(size);}

container::Container3D<Distribution3D::type>& Distribution3D::get_container() {return data;}