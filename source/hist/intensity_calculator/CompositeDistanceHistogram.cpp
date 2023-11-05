#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <constants/Constants.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogram::CompositeDistanceHistogram(
    hist::WeightedDistribution1D&& p_aa, 
    hist::WeightedDistribution1D&& p_aw, 
    hist::WeightedDistribution1D&& p_ww, 
    hist::Distribution1D&& p_tot, 
    const Axis& axis
) : ICompositeDistanceHistogram(std::move(p_tot), axis), p_aa(std::move(p_aa)), p_aw(std::move(p_aw)), p_ww(std::move(p_ww)) {
    use_weighted_sinc_table();
}

CompositeDistanceHistogram::CompositeDistanceHistogram(
    hist::Distribution1D&& p_aa, 
    hist::Distribution1D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::Distribution1D&& p_tot, 
    const Axis& axis
) : ICompositeDistanceHistogram(std::move(p_tot), axis), p_aa(std::move(p_aa)), p_aw(std::move(p_aw)), p_ww(std::move(p_ww)) {}

CompositeDistanceHistogram::~CompositeDistanceHistogram() = default;

const Distribution1D& CompositeDistanceHistogram::get_aa_counts() const {
    return p_aa;
}

Distribution1D& CompositeDistanceHistogram::get_aa_counts() {
    return p_aa;
}

const Distribution1D& CompositeDistanceHistogram::get_aw_counts() const {
    return p_aw;
}

Distribution1D& CompositeDistanceHistogram::get_aw_counts() {
    return p_aw;
}

const Distribution1D& CompositeDistanceHistogram::get_ww_counts() const {
    return p_ww;
}

Distribution1D& CompositeDistanceHistogram::get_ww_counts() {
    return p_ww;
}

void CompositeDistanceHistogram::apply_water_scaling_factor(double k) {
    auto& p_tot = get_total_counts();
    for (unsigned int i = 0; i < get_axis().bins; ++i) {p_tot[i] = p_aa.index(i) + 2*k*p_aw.index(i) + k*k*p_ww.index(i);}
}

auto partial_profile = [] (const Distribution1D& p, const std::vector<double>& q_axis, const view_ptr<const table::DebyeTable> sinqd_table) {
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        Iq[q] = std::inner_product(p.begin(), p.end(), sinqd_table->begin(q), 0.0);
        Iq[q] *= std::exp(-q_axis[q]*q_axis[q]);
    }
    return ScatteringProfile(Iq, debye_axis);
};

const ScatteringProfile CompositeDistanceHistogram::get_profile_aa() const {
    return partial_profile(get_aa_counts(), q_axis, get_sinc_table());
}

const ScatteringProfile CompositeDistanceHistogram::get_profile_aw() const {
    return partial_profile(get_aw_counts(), q_axis, get_sinc_table())*2;
}

const ScatteringProfile CompositeDistanceHistogram::get_profile_ww() const {
    return partial_profile(get_ww_counts(), q_axis, get_sinc_table());
}