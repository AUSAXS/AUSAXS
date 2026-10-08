// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/DistanceHistogram.h>

#include <dataset/Dataset.h>
#include <hist/Histogram.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/HistogramSettings.h>
#include <utility/MultiThreading.h>

#include <numeric>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::hist;

DistanceHistogram::DistanceHistogram() = default;
DistanceHistogram::DistanceHistogram(const DistanceHistogram&) = default;
DistanceHistogram::DistanceHistogram(DistanceHistogram&&) noexcept = default;
DistanceHistogram& DistanceHistogram::operator=(DistanceHistogram&&) noexcept = default;
DistanceHistogram& DistanceHistogram::operator=(const DistanceHistogram&) = default;

DistanceHistogram::DistanceHistogram(hist::Distribution1D&& p_tot) : Histogram( // NOLINT - consumed piecewise below
    std::move(p_tot.get_data()), 
    Axis(0, p_tot.size()*settings::axes::bin_width, p_tot.size())
) {
    initialize();
}

DistanceHistogram::DistanceHistogram(hist::WeightedDistribution1D&& p_tot) : Histogram( // NOLINT - consumed piecewise below
    p_tot.get_content(), 
    Axis(0, p_tot.size()*settings::axes::bin_width, p_tot.size())
) {
    initialize(p_tot.get_weighted_axis());
    sinc_table.set_d_axis(d_axis);
}

DistanceHistogram::DistanceHistogram(std::unique_ptr<ICompositeDistanceHistogram> cdh) : Histogram(cdh->get_counts(), cdh->get_axis()) {
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

template<bool form_factor>
std::vector<double> DistanceHistogram::debye_sum(
    std::span<const double> counts, observer_ptr<const table::DebyeTable> sinqd_table, std::span<const double> q, int first_row
) {
    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(q.size(), 0);
    auto* pool = utility::multi_threading::get_global_pool();
    pool->detach_blocks(0, static_cast<int>(q.size()), // iterate through all q values
        [counts, &Iq, q, first_row, sinqd_table] (int start, int end) {
            for (int i = start; i < end; ++i) {
                Iq[i] = std::transform_reduce(counts.begin(), counts.end(), sinqd_table->begin(first_row+i), 0.0);
                if constexpr (form_factor) {Iq[i] *= std::exp(-q[i]*q[i]);}
            }
        }
    );
    pool->wait();
    return Iq;
}

template<bool form_factor>
ScatteringProfile DistanceHistogram::debye_sum(std::span<const double> counts, observer_ptr<const table::DebyeTable> sinqd_table) {
    Axis debye_axis = constants::axes::q_axis.sub_axis_covering(settings::axes::qmin, settings::axes::qmax);
    int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin
    auto q = std::span<const double>(constants::axes::q_vals).subspan(q0, debye_axis.bins);
    return {debye_sum<form_factor>(counts, sinqd_table, q, q0), debye_axis};
}

ScatteringProfile DistanceHistogram::debye_transform() const {
    return debye_sum<true>(p, sinc_table.get_sinc_table());
}

Dataset DistanceHistogram::debye_transform(const std::vector<double>& q) const {
    return debye_transform<true>(q);
}

template<bool form_factor>
ScatteringProfile DistanceHistogram::debye_transform() const {
    if constexpr (form_factor) {return debye_transform();} // dispatch to any form factors of a subclass
    else {return debye_sum<false>(p, sinc_table.get_sinc_table());}
}

template<bool form_factor>
Dataset DistanceHistogram::debye_transform(const std::vector<double>& q) const {
    // if the q values are within the evaluated default range, we can just interpolate them for better performance
    Axis debye_axis = constants::axes::q_axis.sub_axis_covering(settings::axes::qmin, settings::axes::qmax);
    if (debye_axis.front() <= q.front() && q.back() <= debye_axis.back()) {
        return debye_transform<form_factor>().as_dataset().interpolate(q);
    }
    static table::DebyeTableManager sinc_table_extended;
    sinc_table_extended.set_q_axis(q);
    sinc_table_extended.set_d_axis(this->d_axis);
    return {q, debye_sum<form_factor>(p, sinc_table_extended.get_sinc_table(), q)};
}

template std::vector<double> DistanceHistogram::debye_sum<true>(std::span<const double>, observer_ptr<const table::DebyeTable>, std::span<const double>, int);
template std::vector<double> DistanceHistogram::debye_sum<false>(std::span<const double>, observer_ptr<const table::DebyeTable>, std::span<const double>, int);
template ScatteringProfile DistanceHistogram::debye_sum<true>(std::span<const double>, observer_ptr<const table::DebyeTable>);
template ScatteringProfile DistanceHistogram::debye_sum<false>(std::span<const double>, observer_ptr<const table::DebyeTable>);
template ScatteringProfile DistanceHistogram::debye_transform<true>() const;
template ScatteringProfile DistanceHistogram::debye_transform<false>() const;
template Dataset DistanceHistogram::debye_transform<true>(const std::vector<double>&) const;
template Dataset DistanceHistogram::debye_transform<false>(const std::vector<double>&) const;

const std::vector<double>& DistanceHistogram::get_d_axis() const {return d_axis;}

const std::vector<double>& DistanceHistogram::get_q_axis() {
    static std::vector<double> q_vals; 
    q_vals = constants::axes::q_axis.sub_axis_covering(settings::axes::qmin, settings::axes::qmax).as_vector();
    return q_vals;
}

const std::vector<double>& DistanceHistogram::get_weighted_counts() const {return get_counts();}

bool DistanceHistogram::is_highly_ordered() const {
    return is_highly_ordered(p);
}

bool DistanceHistogram::is_highly_ordered(const std::vector<double>& counts) {
    if (counts.size() < 3) {return false;}

    int peaks = 0;
    int non_zero = 0;
    for (std::size_t i = 1; i + 1 < counts.size(); ++i) {
        if (counts[i] == 0) {continue;}
        if (counts[i] > 1.5*counts[i-1] && counts[i] > 1.5*counts[i+1]) {++peaks;}
        ++non_zero;
    }

    return non_zero != 0 && peaks > non_zero*0.25;
}