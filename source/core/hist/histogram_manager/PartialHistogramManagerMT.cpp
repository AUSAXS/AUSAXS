// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/PartialHistogramManagerMT.h>

#include <data/Molecule.h>
#include <data/state/StateManager.h>  // IWYU pragma: keep
#include <form_factor/FormFactorType.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/CompactCoordinatesFactoryFF.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/histogram_manager/detail/ManagerResults.h>
#include <hist/histogram_manager/detail/PartialBinEstimate.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/HistogramSettings.h>
#include <utility/Logging.h>
#include <utility/MultiThreading.h>


using namespace ausaxs;
using namespace ausaxs::hist;

template<bool weighted_bins, bool form_factors>
PartialHistogramManagerMTBase<weighted_bins, form_factors>::PartialHistogramManagerMTBase(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method) 
    : IPartialHistogramManager(protein), 
      protein(protein),
      exv_method(exv_method),
      coords_a(this->body_size)
{logging::log("initializing PartialHistogramManagerMT");}

template<bool weighted_bins, bool form_factors>
PartialHistogramManagerMTBase<weighted_bins, form_factors>::~PartialHistogramManagerMTBase() = default;

template<bool weighted_bins, bool form_factors>
std::unique_ptr<DistanceHistogram> PartialHistogramManagerMTBase<weighted_bins, form_factors>::calculate() {
    if (!this->statemanager->is_modified() && !cached_p_tot.empty()) {
        logging::log("PartialHistogramManagerMT::calculate: returning cached value");
        auto p_tot = cached_p_tot; // if the state was not modified, we can return the cached value
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
    distance_calculator::Calculator<weighted_bins, form_factors> calculator(*store);

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
    for (int id : recalculated) {store->visit(id, [this] (const auto& result) {hist::detail::fold_classes(this->master, result, std::plus<>());});}
    recalculated.clear();
    this->statemanager->reset_to_false();

    // downsize our axes to only the relevant area
    GenericDistribution1D_t p_tot = this->master; // NOLINT - intentional slicing
    p_tot.resize(hist::detail::trimmed_bin_count(p_tot));

    // update cache
    cached_p_tot = p_tot;

    return std::make_unique<DistanceHistogram>(std::move(p_tot));
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::update_compact_representation_body(int index) {
    const auto& atoms = this->protein->get_body(index).get_atoms();
    if constexpr (form_factors) {
        this->coords_a[index] = hist::detail::factory::construct_by_ff(atoms);
    } else {
        this->coords_a[index] = hist::detail::factory::construct(atoms);
        hist::detail::SimpleExvModel::apply_simple_excluded_volume(this->coords_a[index], this->protein);
    }
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::update_compact_representation_water() {
    this->coords_w = hist::detail::factory::construct_from_waters(this->protein);
}

template<bool weighted_bins, bool form_factors>
std::unique_ptr<ICompositeDistanceHistogram> PartialHistogramManagerMTBase<weighted_bins, form_factors>::calculate_all() {
    logging::log("PartialHistogramManagerMT::calculate_all: starting calculation");
    auto total = calculate();
    int bins = total->get_weighted_counts().size();

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions.
    // only the lower triangle of the body pairs is ever calculated
    using Distributions = hist::detail::ManagerDistributions<weighted_bins, form_factors>;
    std::vector<int> aa_ids;
    for (int i = 0; i < this->body_size; ++i) {
        for (int j = 0; j <= i; ++j) {aa_ids.push_back(aa[i][j]);}
    }

    Distributions d;
    d.p_aa = hist::detail::sum_results<typename Distributions::aa_t>(*store, aa_ids);
    d.p_aw = hist::detail::sum_results<typename Distributions::aw_t>(*store, aw);
    d.p_ww = store->get_1d(ww);
    d.p_tot = this->master; // NOLINT - intentional slicing
    d.resize(bins);
    return hist::detail::make_histogram(std::move(d), exv_method, protein);
}

template<bool weighted_bins, bool form_factors>
int PartialHistogramManagerMTBase<weighted_bins, form_factors>::prepare_axis() {
    int required = hist::detail::required_partial_bin_count(*this->protein);
    if (!this->master.empty()) {
        if (required <= this->master.axis.bins) {return this->master.axis.bins;}

        logging::log("PartialHistogramManagerMT::prepare_axis: structure outgrew its axis; rebuilding");
        this->master = hist::detail::MasterHistogram<weighted_bins>();
        this->statemanager->modified_all();
    }
    return hist::detail::grown_partial_bin_count(required);
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::initialize(int bin_count) {
    Axis axis(0, settings::axes::bin_width*bin_count, bin_count);
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram<weighted_bins>(p_base, axis);
    store = std::make_unique<distance_calculator::HistogramStore<weighted_bins>>(axis.bins, form_factors ? form_factor::get_active_count() : 1);
    aa.assign(this->body_size, std::vector<int>(this->body_size));
    aw.resize(this->body_size);
    for (int n = 0; n < this->body_size; ++n) {
        for (int m = 0; m < this->body_size; ++m) {aa[n][m] = form_factors ? store->allocate_3d() : store->allocate_1d();}
        aw[n] = form_factors ? store->allocate_2d() : store->allocate_1d();
    }
    ww = store->allocate_1d();
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::recalculate(int id) {
    // the result is not written until the calculator runs, so its old contents are still there to be taken out
    store->visit(id, [this] (const auto& result) {hist::detail::fold_classes(this->master, result, std::minus<>());});
    recalculated.push_back(id);
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::calc_self_correlation(calculator_t calculator, int index) {
    update_compact_representation_body(index);
    recalculate(aa[index][index]);
    calculator->enqueue_calculate_self(this->coords_a[index], aa[index][index]);
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::calc_aa(calculator_t calculator, int n, int m) {
    recalculate(aa[n][m]);
    calculator->enqueue_calculate_cross(this->coords_a[n], this->coords_a[m], aa[n][m], 2);
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::calc_aw(calculator_t calculator, int index) {
    recalculate(aw[index]);
    calculator->enqueue_calculate_cross(this->coords_a[index], this->coords_w, aw[index], form_factors ? 1 : 2); // see HistogramManagerMTBase
}

template<bool weighted_bins, bool form_factors>
void PartialHistogramManagerMTBase<weighted_bins, form_factors>::calc_ww(calculator_t calculator) {
    recalculate(ww);
    calculator->enqueue_calculate_self(this->coords_w, ww);
}

template class hist::PartialHistogramManagerMTBase<false, false>;
template class hist::PartialHistogramManagerMTBase<false, true>;
template class hist::PartialHistogramManagerMTBase<true, false>;
template class hist::PartialHistogramManagerMTBase<true, true>;