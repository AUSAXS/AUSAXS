#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <constants/Constants.h>
#include <settings/HistogramSettings.h>

using namespace hist;

template<bool use_weighted_distribution>
CompositeDistanceHistogram<use_weighted_distribution>::CompositeDistanceHistogram(
    typename hist::GenericDistribution1D<use_weighted_distribution>::type&& p_aa, 
    typename hist::GenericDistribution1D<use_weighted_distribution>::type&& p_aw, 
    typename hist::GenericDistribution1D<use_weighted_distribution>::type&& p_ww, 
    hist::Distribution1D&& p_tot, 
    const Axis& axis) 
    : ICompositeDistanceHistogram(std::move(p_tot), axis), p_aa(std::move(p_aa)), p_aw(std::move(p_aw)), p_ww(std::move(p_ww)) 
{}

template<bool use_weighted_distribution>
CompositeDistanceHistogram<use_weighted_distribution>::~CompositeDistanceHistogram() = default;

template<bool use_weighted_distribution>
const std::vector<constants::axes::d_type>& CompositeDistanceHistogram<use_weighted_distribution>::get_aa_counts() const {
    return p_aa.get_counts();
}

template<bool use_weighted_distribution>
std::vector<constants::axes::d_type>& CompositeDistanceHistogram<use_weighted_distribution>::get_aa_counts() {
    return p_aa.get_counts();
}

template<bool use_weighted_distribution>
const std::vector<constants::axes::d_type>& CompositeDistanceHistogram<use_weighted_distribution>::get_aw_counts() const {
    return p_aw.get_counts();
}

template<bool use_weighted_distribution>
std::vector<constants::axes::d_type>& CompositeDistanceHistogram<use_weighted_distribution>::get_aw_counts() {
    return p_aw.get_counts();
}

template<bool use_weighted_distribution>
const std::vector<constants::axes::d_type>& CompositeDistanceHistogram<use_weighted_distribution>::get_ww_counts() const {
    return p_ww.get_counts();
}

template<bool use_weighted_distribution>
std::vector<constants::axes::d_type>& CompositeDistanceHistogram<use_weighted_distribution>::get_ww_counts() {
    return p_ww.get_counts();
}

template<bool use_weighted_distribution>
void CompositeDistanceHistogram<use_weighted_distribution>::apply_water_scaling_factor(double k) {
    auto& p_tot = get_total_counts();
    for (unsigned int i = 0; i < get_axis().bins; ++i) {p_tot[i] = p_aa.index(i) + 2*k*p_aw.index(i) + k*k*p_ww.index(i);}
}

template<bool use_weighted_distribution>
void CompositeDistanceHistogram<use_weighted_distribution>::reset_water_scaling_factor() {
    apply_water_scaling_factor(1);
}

template<bool use_weighted_distribution>
auto partial_profile = [] (const std::vector<constants::axes::d_type>& p, const std::vector<double>& q_axis) {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        Iq[q] = std::inner_product(p.begin(), p.end(), sinqd_table.begin(q), 0.0);
        Iq[q] *= std::exp(-q_axis[q]*q_axis[q]);
    }
    return ScatteringProfile(Iq, debye_axis);
};

template<bool use_weighted_distribution>
const ScatteringProfile CompositeDistanceHistogram<use_weighted_distribution>::get_profile_aa() const {
    return partial_profile<use_weighted_distribution>(get_aa_counts(), q_axis);
}

template<bool use_weighted_distribution>
const ScatteringProfile CompositeDistanceHistogram<use_weighted_distribution>::get_profile_aw() const {
    return partial_profile<use_weighted_distribution>(get_aw_counts(), q_axis)*2;
}

template<bool use_weighted_distribution>
const ScatteringProfile CompositeDistanceHistogram<use_weighted_distribution>::get_profile_ww() const {
    return partial_profile<use_weighted_distribution>(get_ww_counts(), q_axis);
}

template class hist::CompositeDistanceHistogram<false>;
template class hist::CompositeDistanceHistogram<true>;