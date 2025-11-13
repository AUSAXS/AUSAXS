// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <data/state/StateManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/MultiThreading.h>
#include <utility/Logging.h>
#include <container/ThreadLocalWrapper.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool weighted_bins, bool variable_bin_width> 
PartialHistogramManagerMT<weighted_bins, variable_bin_width>::PartialHistogramManagerMT(observer_ptr<const data::Molecule> protein) 
    : PartialHistogramManager<weighted_bins, variable_bin_width>(protein) 
{}

template<bool weighted_bins, bool variable_bin_width> 
PartialHistogramManagerMT<weighted_bins, variable_bin_width>::~PartialHistogramManagerMT() = default;

namespace {
    int water_res_index = 1.31e4;
    int to_res_index(int body1, int body2) {
        return body1 + body2*1e2;
    }

    int to_res_index_water(int body) {
        return body + water_res_index;
    }
}

template<bool weighted_bins, bool variable_bin_width> 
std::unique_ptr<DistanceHistogram> PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calculate() {
    logging::log("PartialHistogramManagerMT::calculate: starting calculation");
    auto& externally_modified = this->statemanager->get_externally_modified_bodies();
    auto& internally_modified = this->statemanager->get_internally_modified_bodies();
    bool hydration_modified = this->statemanager->is_modified_hydration();
    auto pool = utility::multi_threading::get_global_pool();
    auto calculator = std::make_unique<distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>>();

    // check if the object has already been initialized
    if (this->master.empty()) [[unlikely]] {
        initialize(calculator.get()); 

        // since the initialization also calculates the self-correlation, mark it as unmodified to avoid desyncing its state
        internally_modified = std::vector<bool>(this->body_size, false);
    }

    // if not, we must first check if the atom coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < this->body_size; ++i) {

            // if the internal state was modified, we have to recalculate the self-correlation
            if (internally_modified[i]) {
                calc_self_correlation(calculator.get(), i);
            }

            // if the external state was modified, we have to update the coordinate representations for later calculations (implicitly done in calc_self_correlation)
            else if (externally_modified[i]) {
                pool->detach_task(
                    [this, i] () {update_compact_representation_body(i);}
                );
            }
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
        calc_ww(calculator.get());
    }

    // iterate through the lower triangle and check if either of each pair of bodies was modified
    for (unsigned int i = 0; i < this->body_size; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            if (externally_modified[i] || externally_modified[j]) {
                // one of the bodies was modified, so we recalculate its partial histogram
                calc_aa(calculator.get(), i, j);
            }
        }

        // we also have to remember to update the partial histograms with the hydration layer
        if (externally_modified[i] || hydration_modified) {
            calc_aw(calculator.get(), i);
        }
    }

    // merge the partial results from each thread and add it to the master histogram
    // for this process, we first have to wait for all threads to finish
    // then we extract the results in the same order they were submitted to ensure correctness
    auto res = calculator->run();
    {
        if (hydration_modified) {
            assert(res.self.contains(water_res_index) && "PartialHistogramManagerMT::calculate: water result not found");
            pool->detach_task(
                [this, r = std::move(res.self[water_res_index])]
                () mutable {combine_ww(std::move(r));}
            );
        }

        for (unsigned int i = 0; i < this->body_size; ++i) {
            if (internally_modified[i]) {
                assert(res.self.contains(to_res_index(i, i)) && "PartialHistogramManagerMT::calculate: self result not found");
                pool->detach_task(
                    [this, i, r = std::move(res.self[to_res_index(i, i)])]
                    () mutable {combine_self_correlation(i, std::move(r));}
                );
            }

            for (unsigned int j = 0; j < i; ++j) {
                if (externally_modified[i] || externally_modified[j]) {
                    assert(res.cross.contains(to_res_index(i, j)) && "PartialHistogramManagerMT::calculate: cross result not found");
                    pool->detach_task(
                        [this, i, j, r = std::move(res.cross[to_res_index(i, j)])]
                        () mutable {combine_aa(i, j, std::move(r));}
                    );
                }
            }

            if (externally_modified[i] || hydration_modified) {
                assert(res.cross.contains(to_res_index_water(i)) && "PartialHistogramManagerMT::calculate: water result not found");
                pool->detach_task(
                    [this, i, r = std::move(res.cross[to_res_index_water(i)])]
                    () mutable {combine_aw(i, std::move(r));}
                );
            }
        }
    }

    this->statemanager->reset_to_false();
    pool->wait();

    // downsize our axes to only the relevant area
    GenericDistribution1D_t p_tot = this->master;
    int max_bin = 10; // minimum size is 10
    for (int i = (int) p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_tot.resize(max_bin);

    return std::make_unique<DistanceHistogram>(std::move(p_tot));
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::update_compact_representation_body(int index) {
    this->coords_a[index] = detail::CompactCoordinates<variable_bin_width>(this->protein->get_body(index).get_atoms());
    hist::detail::SimpleExvModel::apply_simple_excluded_volume(this->coords_a[index], this->protein);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::update_compact_representation_water() {
    this->coords_w = detail::CompactCoordinates<variable_bin_width>(this->protein->get_waters());
}

template<bool weighted_bins, bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calculate_all() {
    logging::log("PartialHistogramManagerMT::calculate_all: starting calculation");
    auto total = calculate();
    int bins = total->get_total_counts().size();

    // determine p_tot
    GenericDistribution1D_t p_tot(bins);
    for (int i = 0; i < bins; ++i) {
        p_tot.index(i) = this->master.index(i);
    }

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    GenericDistribution1D_t p_ww = this->partials_ww;
    GenericDistribution1D_t p_aa = this->master.base;
    GenericDistribution1D_t p_aw(bins);
    p_ww.resize(bins);
    p_aa.resize(bins);

    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < this->body_size; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            // iterate through each entry in the partial histogram
            std::transform(p_aa.begin(), p_aa.end(), this->partials_aa.index(i, j).begin(), p_aa.begin(), std::plus<>());
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < this->body_size; ++i) {
        // iterate through each entry in the partial histogram
        std::transform(p_aw.begin(), p_aw.end(), this->partials_aw.index(i).begin(), p_aw.begin(), std::plus<>());
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
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::initialize(calculator_t calculator) {
    auto pool = utility::multi_threading::get_global_pool();
    const Axis& axis = constants::axes::d_axis; 
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram<weighted_bins>(p_base, axis);
    this->partials_ww = detail::PartialHistogram<weighted_bins>(axis.bins);
    for (unsigned int i = 0; i < this->body_size; ++i) {
        this->partials_aw.index(i) = detail::PartialHistogram<weighted_bins>(axis.bins);
        this->partials_aa.index(i, i) = detail::PartialHistogram<weighted_bins>(axis.bins);
        calc_self_correlation(calculator, i);

        for (unsigned int j = 0; j < i; ++j) {
            this->partials_aa.index(i, j) = detail::PartialHistogram<weighted_bins>(axis.bins);
        }
    }

    auto res = calculator->run();
    assert(res.self.size() == this->body_size && "The number of self-correlation results does not match the number of bodies.");
    for (int i = 0; i < static_cast<int>(this->body_size); ++i) {
        assert(res.self.contains(to_res_index(i, i)) && "PartialHistogramManagerMT::initialize: self result not found");
        pool->detach_task(
            [this, i, r = std::move(res.self[to_res_index(i, i)])] () mutable {combine_self_correlation(i, std::move(r));}
        );
    }
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_self_correlation(calculator_t calculator, int index) {
    update_compact_representation_body(index);
    calculator->enqueue_calculate_self(this->coords_a[index], 1, to_res_index(index, index));
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_aa(calculator_t calculator, int n, int m) {
    calculator->enqueue_calculate_cross(this->coords_a[n], this->coords_a[m], 1, to_res_index(n, m));
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_aw(calculator_t calculator, int index) {
    calculator->enqueue_calculate_cross(this->coords_a[index], this->coords_w, 1, to_res_index_water(index));
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::calc_ww(calculator_t calculator) {
    calculator->enqueue_calculate_self(this->coords_w, 1, water_res_index);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::combine_self_correlation(int index, GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(index, index);
    this->partials_aa.index(index, index) = std::move(res);
    this->master += this->partials_aa.index(index, index);
    master_hist_mutex.unlock();
}

template<bool weighted_bins, bool variable_bin_width> 
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::combine_aa(int n, int m, GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(n, m);
    this->partials_aa.index(n, m) = std::move(res);
    this->master += this->partials_aa.index(n, m);
    master_hist_mutex.unlock();
}

template<bool weighted_bins, bool variable_bin_width> 
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::combine_aw(int index, GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_aw.index(index);
    this->partials_aw.index(index) = std::move(res);
    this->master += this->partials_aw.index(index);
    master_hist_mutex.unlock();
}

template<bool weighted_bins, bool variable_bin_width> 
void PartialHistogramManagerMT<weighted_bins, variable_bin_width>::combine_ww(GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_ww;
    this->partials_ww = std::move(res);
    this->master += this->partials_ww;
    master_hist_mutex.unlock();
}

template class hist::PartialHistogramManagerMT<false, false>;
template class hist::PartialHistogramManagerMT<false, true>;
template class hist::PartialHistogramManagerMT<true, false>;
template class hist::PartialHistogramManagerMT<true, true>;