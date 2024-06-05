/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <settings/HistogramSettings.h>

using namespace hist;

Limit CompositeDistanceHistogramFoXS::get_excluded_volume_scaling_factor_limits() const {
    return {0.8, 1.2};
}

double CompositeDistanceHistogramFoXS::exv_factor(double q) const {
    constexpr double rm = 1.58;
    constexpr double c = rm*rm/(4*constants::pi);
    return std::pow(free_params.cx, 3)*std::exp(-c*(std::pow(free_params.cx, 2) - 1)*q*q);
}

CompositeDistanceHistogramFoXS::CompositeDistanceHistogramFoXS(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot
) : CompositeDistanceHistogramFFExplicitBase(std::move(p_aa), std::move(p_ax), std::move(p_xx), std::move(p_aw), std::move(p_wx), std::move(p_ww), std::move(p_tot)) {
    /* 
     * All p(r) calculators use the charge of the atom as a weight for the histogram, which are then multiplied by normalized form factors here.
     * Only the excluded volume dummy atoms use unity weights in our implementation.
     * 
     * FoXS flips this around, calculating a histogram using unity weights for each interaction type, and then multiplying by their scaled form factors here. 
     * Thus we redefine the weighted histograms with the excluded volume histogram, since that is the only one using unity weights. 
     */
    distance_profiles.aa = exv_distance_profiles.xx;
    distance_profiles.aw = exv_distance_profiles.wx;
}

CompositeDistanceHistogramFoXS::CompositeDistanceHistogramFoXS(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : CompositeDistanceHistogramFFExplicitBase(std::move(p_aa), std::move(p_ax), std::move(p_xx), std::move(p_aw), std::move(p_wx), std::move(p_ww), std::move(p_tot)) {
    distance_profiles.aa = exv_distance_profiles.xx;
    distance_profiles.aw = exv_distance_profiles.wx;
}