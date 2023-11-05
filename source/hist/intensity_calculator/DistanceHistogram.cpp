#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distribution/WeightedDistribution.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <table/VectorDebyeTable.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>
#include <constants/Constants.h>

using namespace hist;

DistanceHistogram::DistanceHistogram() = default;

DistanceHistogram::DistanceHistogram(hist::Distribution1D&& p_tot, const Axis& axis) : Histogram(std::move(p_tot.get_data()), axis) {
    initialize();
}

DistanceHistogram::DistanceHistogram(hist::WeightedDistribution1D&& p_tot, const Axis& axis) : Histogram(std::move(p_tot.get_data()), axis) {
    use_weighted_sinc_table();
    initialize();
}

DistanceHistogram::DistanceHistogram(std::unique_ptr<ICompositeDistanceHistogram> cdh) : Histogram(std::move(cdh->get_counts()), cdh->get_axis()) {
    initialize();
}

const view_ptr<const table::DebyeTable> DistanceHistogram::get_sinc_table() const {
    if (use_weighted_table) {return weighted_sinc_table;}
    return view_ptr<const table::DebyeTable>(table::ArrayDebyeTable::get_default_table());
}

void DistanceHistogram::use_weighted_sinc_table() {
    auto weighted_bins = hist::WeightedDistribution::get_weighted_bins();
    weighted_sinc_table = std::make_unique<table::VectorDebyeTable>(table::VectorDebyeTable(weighted_bins));
    use_weighted_table = true;
}

DistanceHistogram::~DistanceHistogram() = default;

void DistanceHistogram::initialize() {
    d_axis = axis.as_vector();
    d_axis[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms
    table::ArrayDebyeTable::check_default(d_axis);
}

ScatteringProfile DistanceHistogram::debye_transform() const {
    // calculate the Debye scattering intensity
    const auto& q_axis = constants::axes::q_vals;
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    const auto& sinqd_table = get_sinc_table();

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) { // iterate through all q values
        Iq[q] = std::inner_product(p.begin(), p.end(), sinqd_table->begin(q), 0.0);
        Iq[q] *= std::exp(-q_axis[q]*q_axis[q]); // form factor
    }
    return ScatteringProfile(Iq, debye_axis);
}

const std::vector<double>& DistanceHistogram::get_d_axis() const {return d_axis;}

const std::vector<double>& DistanceHistogram::get_q_axis() {
    static const std::vector<double> q_vals = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
    return q_vals;
}

const std::vector<double>& DistanceHistogram::get_total_counts() const {return get_counts();}

std::vector<double>& DistanceHistogram::get_total_counts() {return get_counts();}