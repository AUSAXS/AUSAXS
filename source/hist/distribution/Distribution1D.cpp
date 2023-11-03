#include <hist/distribution/Distribution1D.h>

#include <cmath>

using namespace hist;

void Distribution1D::add(float distance, constants::axes::d_type value) {index(std::round(distance)) += value;}
void Distribution1D::add(int32_t i, constants::axes::d_type value) {index(i) += value;}