/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>

using namespace hist;

Limit CompositeDistanceHistogramFoXS::get_excluded_volume_scaling_factor_limits() const {
    return {0.95, 1.05};
}

double CompositeDistanceHistogramFoXS::exv_factor(double q) const {
    constexpr double rm = 1.58;
    constexpr double c = rm*rm/(4*constants::pi);
    return std::pow(cx, 3)*std::exp(-c*(cx*cx - 1)*q*q);
}