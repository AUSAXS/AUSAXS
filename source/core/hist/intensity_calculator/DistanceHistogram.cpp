// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <table/VectorDebyeTable.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace ausaxs::hist;

DistanceHistogram::DistanceHistogram() = default;
DistanceHistogram::DistanceHistogram(const DistanceHistogram&) = default;
DistanceHistogram::DistanceHistogram(DistanceHistogram&&) noexcept = default;
DistanceHistogram& DistanceHistogram::operator=(DistanceHistogram&&) noexcept = default;
DistanceHistogram& DistanceHistogram::operator=(const DistanceHistogram&) = default;

DistanceHistogram::DistanceHistogram(hist::Distribution1D&& p_tot) : Histogram(
    std::move(p_tot.get_data()), 
    Axis(0, p_tot.size()*settings::axes::bin_width, p_tot.size())
) {
    initialize();
}

DistanceHistogram::DistanceHistogram(hist::WeightedDistribution1D&& p_tot) : Histogram(
    p_tot.get_content(), 
    Axis(0, p_tot.size()*settings::axes::bin_width, p_tot.size())
) {
    initialize(p_tot.get_weighted_axis());
    sinc_table.set_d_axis(d_axis);
}

DistanceHistogram::DistanceHistogram(std::unique_ptr<ICompositeDistanceHistogram> cdh) : Histogram(std::move(cdh->get_counts()), cdh->get_axis()) {
    initialize();
}

DistanceHistogram::~DistanceHistogram() = default;

void DistanceHistogram::initialize(std::vector<double>&& d_axis) {
    this->d_axis = std::move(d_axis);
    this->d_axis[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms
}

void DistanceHistogram::initialize() {
    d_axis = axis.as_vector();
    d_axis[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms
}

ScatteringProfile DistanceHistogram::debye_transform() const {
    // calculate the Debye scattering intensity
    const auto& q_axis = constants::axes::q_vals;
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    auto sinqd_table = sinc_table.get_sinc_table();

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin
    for (int q = q0; q < static_cast<int>(q0+debye_axis.bins); ++q) { // iterate through all q values
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
    sinc_table.set_q_axis(q);
    const auto& sinqd_table = sinc_table.get_sinc_table();

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(q.size(), 0);
    for (int i = 0; i < static_cast<int>(q.size()); ++i) { // iterate through all q values
        Iq[i] = std::inner_product(p.begin(), p.end(), sinqd_table->begin(i), 0.0);
        Iq[i] *= std::exp(-q[i]*q[i]); // form factor
    }
    return SimpleDataset(q, Iq);
}

const std::vector<double>& DistanceHistogram::get_d_axis() const {return d_axis;}

const std::vector<double>& DistanceHistogram::get_q_axis() {
    static std::vector<double> q_vals; 
    q_vals = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
    return q_vals;
}

const std::vector<double>& DistanceHistogram::get_weighted_counts() const {return get_counts();}

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