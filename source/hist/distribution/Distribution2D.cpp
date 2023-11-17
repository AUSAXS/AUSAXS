#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution2D.h>

#include <cmath>

using namespace hist;

Distribution2D::Distribution2D(WeightedDistribution2D&& other) : Container2D(std::move(other)) {}
void Distribution2D::add(unsigned int x, float distance, constants::axes::d_type value) {index(x, std::round(distance)) += value;}
void Distribution2D::add(unsigned int x, int32_t i, constants::axes::d_type value) {index(x, i) += value;}