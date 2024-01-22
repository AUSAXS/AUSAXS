#include <hist/distribution/Distribution2D.h>

#include <cmath>

using namespace hist;

void Distribution2D::add(unsigned int x, float distance, constants::axes::d_type value) {index(x, std::round(distance)) += value;}
void Distribution2D::add(unsigned int x, int32_t i, constants::axes::d_type value) {index(x, i) += value;}