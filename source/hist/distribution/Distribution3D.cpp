#include <hist/distribution/Distribution3D.h>

#include <cmath>

using namespace hist;

void Distribution3D::add(unsigned int x, unsigned int y, float distance, constants::axes::d_type value) {index(x, y, std::round(distance)) += value;}
void Distribution3D::add(unsigned int x, unsigned int y, int32_t i, constants::axes::d_type value) {index(x, y, i) += value;}