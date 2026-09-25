// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/PartialHistogramManagerMT.h>

#include <data/Molecule.h>
#include <data/state/StateManager.h>  // IWYU pragma: keep
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <settings/HistogramSettings.h>
#include <utility/Logging.h>
#include <utility/MultiThreading.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool weighted_bins, bool variable_bin_width> 
PartialHistogramManagerMT<weighted_bins, variable_bin_width>::PartialHistogramManagerMT(observer_ptr<const data::Molecule> protein) 
    : PartialHistogramManager<weighted_bins, variable_bin_width>(protein) 
{logging::log("initializing PartialHistogramManagerMT");}

template<bool weighted_bins, bool variable_bin_width> 
PartialHistogramManagerMT<weighted_bins, variable_bin_width>::~PartialHistogramManagerMT() = default;

template<bool weighted_bins, bool variable_bin_width> 
std::unique_ptr<DistanceHistogram> PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calculate() {
    if (!this->statemanager->is_modified() && !cache.p_tot.empty()) {
        logging::log("PartialHistogramManagerMT::calculate: returning cached value");
        auto p_tot = cache.p_tot; // if the state was not modified, we can return the cached value
        return std::make_unique<DistanceHistogram>(std::move(p_tot));
    }

    logging::log("PartialHistogramManagerMT::calculate: starting calculation");
    int bin_count = this->prepare_axis();
    const auto& externally_modified = this->statemanager->get_externally_modified_bodies();
    const auto& internally_modified = this->statemanager->get_internally_modified_bodies();
    bool hydration_modified = this->statemanager->is_modified_hydration();
    auto* pool = utility::multi_threading::get_global_pool();

    // check if the object has already been initialized
    bool initialized = !this->master.empty();
    if (!initialized) [[unlikely]] {
        initialize(bin_count);
    }
    distance_calculator::Calculator<weighted_bins, variable_bin_width> calculator(*store);

    for (int i = 0; i < this->body_size; ++i) {
        // the self-correlation is calculated when the body is first seen, and again whenever its internal state was modified
        if (!initialized || internally_modified[i]) {
            calc_self_correlation(&calculator, i);
        }

        // if only the external state was modified, we have to update the coordinate representations for later calculations (implicitly done in calc_self_correlation)
        else if (externally_modified[i]) {
            pool->detach_task(
                [this, i] () {update_compact_representation_body(i);}
            );
        }
    }

    // small efficiency improvement: if the hydration layer was modified, we can update the compact representations in parallel with the self-correlation
    if (hydration_modified) {
        pool->detach_task(
            [this] () {update_compact_representation_water();}
        );
    }
    pool->wait(); // ensure the compact representations have been updated before continuing

    // check if the hydration layer was modified
    if (hydration_modified) {
        calc_ww(&calculator);
    }

    // iterate through the lower triangle and check if either of each pair of bodies was modified
    for (int i = 0; i < this->body_size; ++i) {
        // everything body i pairs with is enqueued back-to-back with no work in between, so it is held
        // and dispatched as one group; the coordinates were all built above, before pool->wait()
        calculator.hold();
        for (int j = 0; j < i; ++j) {
            if (externally_modified[i] || externally_modified[j]) {
                // one of the bodies was modified, so we recalculate its partial histogram
                calc_aa(&calculator, i, j);
            }
        }

        // we also have to remember to update the partial histograms with the hydration layer
        if (externally_modified[i] || hydration_modified) {
            calc_aw(&calculator, i);
        }
        calculator.release_hold();
    }

    // the recalculated partial histograms replace their old contents in the store, which were taken out of the master histogram as they were queued
    calculator.run();
    for (int h : recalculated) {this->master += store->row(h);}
    recalculated.clear();
    this->statemanager->reset_to_false();

    // downsize our axes to only the relevant area
    GenericDistribution1D_t p_tot = this->master; // NOLINT - intentional slicing
    p_tot.resize(hist::detail::trimmed_bin_count(p_tot));

    // update cache
    cache.p_tot = p_tot;

    return std::make_unique<DistanceHistogram>(std::move(p_tot));
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::update_compact_representation_body(int index) {
    this->coords_a[index] = hist::detail::factory::construct<variable_bin_width>(this->protein->get_body(index).get_atoms());
    hist::detail::SimpleExvModel::apply_simple_excluded_volume(this->coords_a[index], this->protein);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::update_compact_representation_water() {
    this->coords_w = hist::detail::factory::construct_from_waters<variable_bin_width>(this->protein);
}

template<bool weighted_bins, bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calculate_all() {
    if (
        !this->statemanager->is_modified() 
        && !cache.p_tot.empty() && !cache.p_aa.empty() && !cache.p_aw.empty() && !cache.p_ww.empty()
    ) {
        logging::log("PartialHistogramManagerMT::calculate_all: returning cached value");
        auto p_tot = cache.p_tot; // if the state was not modified, we can return the cached value
        auto p_aa = cache.p_aa;
        auto p_aw = cache.p_aw;
        auto p_ww = cache.p_ww;
        if constexpr (weighted_bins) {
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

    logging::log("PartialHistogramManagerMT::calculate_all: starting calculation");
    auto total = calculate();
    int bins = total->get_weighted_counts().size();

    // determine p_tot
    GenericDistribution1D_t p_tot(bins);
    for (int i = 0; i < bins; ++i) {
        p_tot.index(i) = this->master.index(i);
    }

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    GenericDistribution1D_t p_ww = store->export_1d(handle_ww());
    GenericDistribution1D_t p_aa = this->master.base;
    GenericDistribution1D_t p_aw(bins);
    p_ww.resize(bins);
    p_aa.resize(bins);

    // iterate through all partial histograms in the upper triangle
    for (int i = 0; i < this->body_size; ++i) {
        for (int j = 0; j <= i; ++j) {
            // iterate through each entry in the partial histogram
            std::transform(p_aa.begin(), p_aa.end(), store->row(handle_aa(i, j)).begin(), p_aa.begin(), std::plus<>());
        }
    }

    // iterate through all partial hydration-protein histograms
    for (int i = 0; i < this->body_size; ++i) {
        // iterate through each entry in the partial histogram
        std::transform(p_aw.begin(), p_aw.end(), store->row(handle_aw(i)).begin(), p_aw.begin(), std::plus<>());
    }

    if constexpr (weighted_bins) {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(Distribution1D(p_aa)), 
            std::move(Distribution1D(p_aw)), 
            std::move(Distribution1D(p_ww)), 
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

template<bool weighted_bins, bool variable_bin_width>
int PartialHistogramManagerMT<weighted_bins, variable_bin_width>::handle_aa(int n, int m) const {
    assert(0 <= m && m <= n && n < this->body_size && "PartialHistogramManagerMT::handle_aa: expected a body pair in the lower triangle.");
    return n*(n+1)/2 + m;
}

template<bool weighted_bins, bool variable_bin_width>
int PartialHistogramManagerMT<weighted_bins, variable_bin_width>::handle_aw(int index) const {
    return this->body_size*(this->body_size+1)/2 + index;
}

template<bool weighted_bins, bool variable_bin_width>
int PartialHistogramManagerMT<weighted_bins, variable_bin_width>::handle_ww() const {
    return handle_aw(this->body_size);
}

template<bool weighted_bins, bool variable_bin_width> 
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::initialize(int bin_count) {
    Axis axis(0, settings::axes::bin_width*bin_count, bin_count);
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram<weighted_bins>(p_base, axis);
    store = std::make_unique<distance_calculator::HistogramStore<weighted_bins>>(handle_ww()+1, axis.bins);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::recalculate(int h) {
    // the row is not written until the calculator runs, so its old contents are still there to be taken out
    this->master -= store->row(h);
    recalculated.push_back(h);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_self_correlation(calculator_t calculator, int index) {
    update_compact_representation_body(index);
    recalculate(handle_aa(index, index));
    calculator->enqueue_calculate_self(this->coords_a[index], handle_aa(index, index));
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_aa(calculator_t calculator, int n, int m) {
    recalculate(handle_aa(n, m));
    calculator->enqueue_calculate_cross(this->coords_a[n], this->coords_a[m], handle_aa(n, m), 2);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_aw(calculator_t calculator, int index) {
    recalculate(handle_aw(index));
    calculator->enqueue_calculate_cross(this->coords_a[index], this->coords_w, handle_aw(index), 2);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_ww(calculator_t calculator) {
    recalculate(handle_ww());
    calculator->enqueue_calculate_self(this->coords_w, handle_ww());
}

template class hist::PartialHistogramManagerMT<false, false>;
template class hist::PartialHistogramManagerMT<false, true>;
template class hist::PartialHistogramManagerMT<true, false>;
template class hist::PartialHistogramManagerMT<true, true>;