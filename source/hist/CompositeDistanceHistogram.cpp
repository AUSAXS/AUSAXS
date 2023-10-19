#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>

using namespace hist;

CompositeDistanceHistogram::CompositeDistanceHistogram(std::vector<double>&& p_pp, std::vector<double>&& p_hp, std::vector<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis) 
    : DistanceHistogram(std::move(p_tot), axis), p_aa(std::move(p_aa)), p_wa(std::move(p_wa)), p_ww(std::move(p_ww)) {}

CompositeDistanceHistogram::CompositeDistanceHistogram(std::vector<double>&& p_tot, const Axis& axis) : DistanceHistogram(std::move(p_tot), axis) {}

CompositeDistanceHistogram::~CompositeDistanceHistogram() = default;

const std::vector<double>& CompositeDistanceHistogram::get_pp_counts() const {return p_aa;}

const std::vector<double>& CompositeDistanceHistogram::get_hp_counts() const {return p_wa;}

const std::vector<double>& CompositeDistanceHistogram::get_hh_counts() const {return p_ww;}

void CompositeDistanceHistogram::apply_water_scaling_factor(double k) {
    auto& p_tot = get_total_counts();
    for (unsigned int i = 0; i < get_axis().bins; ++i) {p_tot[i] = p_aa[i] + 2*k*p_wa[i] + k*k*p_ww[i];}
}

void CompositeDistanceHistogram::reset_water_scaling_factor() {
    apply_water_scaling_factor(1);
}