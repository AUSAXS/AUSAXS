// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <data/Molecule.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool use_weighted_distribution>
HistogramManagerMT<use_weighted_distribution>::~HistogramManagerMT() = default;

template<bool use_weighted_distribution>
std::unique_ptr<DistanceHistogram> HistogramManagerMT<use_weighted_distribution>::calculate() {return calculate_all();}

template<bool use_weighted_distribution>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMT<use_weighted_distribution>::calculate_all() {
    logging::log("HistogramManagerMT::calculate: starting calculation");
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time (class is not meant for serial use anyway?)
    hist::detail::CompactCoordinates data_a(this->protein->get_bodies());
    hist::detail::CompactCoordinates data_w(this->protein->get_waters());
    hist::detail::SimpleExvModel::apply_simple_excluded_volume(data_a, this->protein);

    auto calculator = distance_calculator::SimpleCalculator<use_weighted_distribution>::create();
    calculator->enqueue_calculate_self(data_a);
    calculator->enqueue_calculate_self(data_w);
    calculator->enqueue_calculate_cross(data_a, data_w);
    auto res = calculator->run();

    auto p_aa = res.self[0];
    auto p_ww = res.self[1];
    auto p_aw = res.cross[0];

    // calculate p_tot
    GenericDistribution1D_t p_tot(constants::axes::d_axis.bins);
    for (unsigned int i = 0; i < p_tot.size(); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + p_aw.index(i);}

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_aa.resize(max_bin);
    p_ww.resize(max_bin);
    p_aw.resize(max_bin);
    p_tot.resize(max_bin);

    if constexpr (use_weighted_distribution) {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(Distribution1D(std::move(p_aa))), 
            std::move(Distribution1D(std::move(p_aw))), 
            std::move(Distribution1D(std::move(p_ww))), 
            std::move(p_tot)
        );
    } else {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(p_aa), 
            std::move(p_aw), 
            std::move(p_ww), 
            std::move(p_tot)
        );
    }
}

template class hist::HistogramManagerMT<false>;
template class hist::HistogramManagerMT<true>;