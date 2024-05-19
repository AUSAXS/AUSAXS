/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "constants/Axes.h"
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <table/VectorDebyeTable.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>
#include <constants/Constants.h>

using namespace hist;

DistanceHistogram::DistanceHistogram() = default;

DistanceHistogram::DistanceHistogram(DistanceHistogram&& other) = default;

DistanceHistogram::DistanceHistogram(hist::Distribution1D&& p_tot) : Histogram(std::move(p_tot.get_data()), Axis(0, p_tot.size()*constants::axes::d_axis.width(), p_tot.size())) {
    initialize();
}

DistanceHistogram::DistanceHistogram(hist::WeightedDistribution1D&& p_tot) : Histogram(p_tot.get_content(), Axis(0, p_tot.size()*constants::axes::d_axis.width(), p_tot.size())) {
    initialize(p_tot.get_weighted_axis());
    use_weighted_sinc_table();
}

DistanceHistogram::DistanceHistogram(std::unique_ptr<ICompositeDistanceHistogram> cdh) : Histogram(std::move(cdh->get_counts()), cdh->get_axis()) {
    initialize();
}

observer_ptr<const table::DebyeTable> DistanceHistogram::get_sinc_table() const {
    if (use_weighted_table) {return weighted_sinc_table.get();}
    return &table::ArrayDebyeTable::get_default_table();
}

void DistanceHistogram::use_weighted_sinc_table() {
    weighted_sinc_table = std::make_unique<table::VectorDebyeTable>(d_axis);
    use_weighted_table = true;
}

DistanceHistogram::~DistanceHistogram() = default;

void DistanceHistogram::initialize(std::vector<double>&& d_axis) {
    this->d_axis = std::move(d_axis);
    this->d_axis[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms
}

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
        Iq[q-q0] = std::inner_product(p.begin(), p.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] *= std::exp(-q_axis[q]*q_axis[q]); // form factor
    }
    return ScatteringProfile(Iq, debye_axis);
}

SimpleDataset DistanceHistogram::debye_transform(const std::vector<double>& q) const {
    // if the q values are within the default range, we can just interpolate them for better performance
    if (constants::axes::q_axis.min < q.front() && q.back() < constants::axes::q_axis.max) {
        return debye_transform().as_dataset().interpolate(q);
    }

    auto sinqd_table = table::VectorDebyeTable(d_axis, q);

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(q.size(), 0);
    for (unsigned int i = 0; i < q.size(); ++i) { // iterate through all q values
        Iq[i] = std::inner_product(p.begin(), p.end(), sinqd_table.begin(i), 0.0);
        Iq[i] *= std::exp(-q[i]*q[i]); // form factor
    }
    return SimpleDataset(q, Iq);
}

const std::vector<double>& DistanceHistogram::get_d_axis() const {return d_axis;}

const std::vector<double>& DistanceHistogram::get_q_axis() {
    static const std::vector<double> q_vals = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
    return q_vals;
}

const std::vector<double>& DistanceHistogram::get_total_counts() const {return get_counts();}

std::vector<double>& DistanceHistogram::get_total_counts() {return get_counts();}

bool DistanceHistogram::is_highly_ordered() const {
    // count the number of 'spikes', defined as points that are at least 50% larger than their neighbours
    // also count the number of bins with non-zero counts
    unsigned int peaks = 0;
    unsigned int non_zero = 0;
    for (unsigned int i = 1; i < p.size()-1; ++i) {
        if (p[i] == 0) {continue;}
        if (p[i] > 1.5*p[i-1] && p[i] > 1.5*p[i+1]) {++peaks;}
        ++non_zero;
    }

    // if the number of peaks is at least 25% of the number of non-zero bins, we consider the structure to be highly ordered
    return peaks > non_zero*0.25;
}