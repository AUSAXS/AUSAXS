#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <constants/Constants.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogram::CompositeDistanceHistogram(std::vector<double>&& p_aa, std::vector<double>&& p_aw, std::vector<double>&& p_ww, std::vector<double>&& p_tot, const Axis& axis) 
    : ICompositeDistanceHistogram(std::move(p_tot), axis), p_aa(std::move(p_aa)), p_aw(std::move(p_aw)), p_ww(std::move(p_ww)) {}

CompositeDistanceHistogram::CompositeDistanceHistogram(std::vector<double>&& p_tot, const Axis& axis) : ICompositeDistanceHistogram(std::move(p_tot), axis) {}

CompositeDistanceHistogram::~CompositeDistanceHistogram() = default;

const std::vector<double>& CompositeDistanceHistogram::get_aa_counts() const {return p_aa;}
std::vector<double>& CompositeDistanceHistogram::get_aa_counts() {return p_aa;}

const std::vector<double>& CompositeDistanceHistogram::get_aw_counts() const {return p_aw;}
std::vector<double>& CompositeDistanceHistogram::get_aw_counts() {return p_aw;}

const std::vector<double>& CompositeDistanceHistogram::get_ww_counts() const {return p_ww;}
std::vector<double>& CompositeDistanceHistogram::get_ww_counts() {return p_ww;}

void CompositeDistanceHistogram::apply_water_scaling_factor(double k) {
    auto& p_tot = get_total_counts();
    for (unsigned int i = 0; i < get_axis().bins; ++i) {p_tot[i] = p_aa[i] + 2*k*p_aw[i] + k*k*p_ww[i];}
}

void CompositeDistanceHistogram::reset_water_scaling_factor() {
    apply_water_scaling_factor(1);
}

auto partial_profile = [] (const std::vector<double>& p, const std::vector<double>& q_axis) {
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

const ScatteringProfile CompositeDistanceHistogram::get_profile_aa() const {
    return partial_profile(p_aa, q_axis);
}

const ScatteringProfile CompositeDistanceHistogram::get_profile_aw() const {
    return partial_profile(p_aw, q_axis)*2;
}

const ScatteringProfile CompositeDistanceHistogram::get_profile_ww() const {
    return partial_profile(p_ww, q_axis);
}