#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>

using namespace hist;

CompositeDistanceHistogram::CompositeDistanceHistogram(std::vector<double>&& p_pp, std::vector<double>&& p_hp, std::vector<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis) 
    : DistanceHistogram(std::move(p_tot), axis), p_pp(std::move(p_pp)), p_hp(std::move(p_hp)), p_hh(std::move(p_hh)) {}

CompositeDistanceHistogram::CompositeDistanceHistogram(std::vector<double>&& p_tot, const Axis& axis) : DistanceHistogram(std::move(p_tot), axis) {}

CompositeDistanceHistogram::~CompositeDistanceHistogram() = default;

const std::vector<double>& CompositeDistanceHistogram::get_pp_counts() const {return p_pp;}

const std::vector<double>& CompositeDistanceHistogram::get_hh_counts() const {return p_hh;}

const std::vector<double>& CompositeDistanceHistogram::get_hp_counts() const {return p_hp;}

void CompositeDistanceHistogram::apply_water_scaling_factor(double k) {
    auto& p_tot = get_total_counts();
    for (unsigned int i = 0; i < get_axis().bins; ++i) {p_tot[i] = p_pp[i] + 2*k*p_hp[i] + k*k*p_hh[i];}
}

void CompositeDistanceHistogram::reset_water_scaling_factor() {
    apply_water_scaling_factor(1);
}