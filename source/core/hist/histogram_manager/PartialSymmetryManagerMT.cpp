// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/PartialSymmetryManagerMT.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/state/StateManager.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <form_factor/FormFactorType.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/histogram_manager/detail/ManagerResults.h>
#include <hist/histogram_manager/detail/PartialBinEstimate.h>
#include <hist/histogram_manager/detail/SymmetryHelpers.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/HistogramSettings.h>
#include <utility/Logging.h>
#include <utility/MultiThreading.h>

#include <cassert>

/**
The indexing in this file is a bit tricky. 
The body symmetries (body.symmetry.get) contains the first actual symmetry at index 0
The coordinates are indexed such that the main body is at index 0, with the symmetries starting from index 1
Thus, we have a lot of +1s and -1s in the indexing to account for this.
**/

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool weighted_bins, bool form_factors>
PartialSymmetryManagerMTBase<weighted_bins, form_factors>::PartialSymmetryManagerMTBase(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method) 
    : IPartialHistogramManager(protein), 
      protein(protein),
      exv_method(exv_method),
      coords(this->body_size)
{}

template<bool weighted_bins, bool form_factors>
PartialSymmetryManagerMTBase<weighted_bins, form_factors>::~PartialSymmetryManagerMTBase() = default;

template<bool weighted_bins, bool form_factors>
std::unique_ptr<DistanceHistogram> PartialSymmetryManagerMTBase<weighted_bins, form_factors>::calculate() {
    logging::log("PartialSymmetryManagerMT::calculate: starting calculation");
    if (protein->size_water() == 0 && !this->statemanager->is_modified_hydration()) {
        return _calculate<false>();
    }
    return _calculate<true>();
}

template<bool weighted_bins, bool form_factors> template<bool hydration_enabled>
std::unique_ptr<DistanceHistogram> PartialSymmetryManagerMTBase<weighted_bins, form_factors>::_calculate() {
    if (!this->statemanager->is_modified() && !cached_p_tot.empty()) {
        logging::log("PartialSymmetryManagerMT::calculate: returning cached value");
        auto p_tot = cached_p_tot; // if the state was not modified, we can return the cached value
        return std::make_unique<DistanceHistogram>(std::move(p_tot));
    }

    int bin_count = prepare_axis();
    auto externally_modified = this->statemanager->get_externally_modified_bodies();
    auto internally_modified = this->statemanager->get_internally_modified_bodies();
    auto symmetry_modified = this->statemanager->get_symmetry_modified_bodies();
    bool hydration_modified = this->statemanager->is_modified_hydration();

    // shared reference symmetries couple several bodies; expand the flags so that a change to any
    // participating body (or the shared symmetry) recomputes the whole group's copies
    propagate_reference_symmetry_modifications(externally_modified, internally_modified, symmetry_modified);

    auto* pool = utility::multi_threading::get_global_pool();

    // check if the object has already been initialized
    if (this->master.empty()) [[unlikely]] {
        initialize(bin_count);
    }

    // if not, we must first check if the atom coordinates have been changed in any of the bodies
    else {
        for (int ibody = 0; ibody < this->body_size; ++ibody) {

            // if the internal state was modified, we have to recalculate the self-correlation
            if (internally_modified[ibody]) {
                update_compact_representation_body(ibody); //? unnecessary to update whole body; enough to update main body
            }

            // if the external state was modified, we have to update the coordinate representations for later calculations
            // (implicitly done in calc_self_correlation)
            else if (externally_modified[ibody]) {
                pool->detach_task(
                    [this, ibody] () {update_compact_representation_body(ibody);}
                );
            }

            // only update individual symmetry copies when the body itself is not being fully regenerated;
            // concurrent body + symmetry updates on the same body would race on coords[ibody]
            else {
                for (int isym = 0; isym < this->protein->get_body(ibody).size_symmetry(); ++isym) {
                    if (symmetry_modified[ibody][isym]) {
                        pool->detach_task(
                            [this, ibody, isym] () {update_compact_representation_symmetry(ibody, isym+1);}
                        );
                    }
                }
            }
        }
    }

    if constexpr (hydration_enabled) {
        // small efficiency improvement: if the hydration layer was modified, 
        // we can update the compact representations in parallel with the self-correlation
        if (hydration_modified) {
            pool->detach_task(
                [this] () {update_compact_representation_water();}
            );
        }
    }
    pool->wait(); // ensure the compact representations have been updated before continuing

    distance_calculator::Calculator<weighted_bins, form_factors> calculator(*store);

    if constexpr (hydration_enabled) {
        // check if the hydration layer was modified
        if (hydration_modified) {
            calc_ww(&calculator);
        }
    }

    // internal modification implies external modification, which in turn moves every copy of the body
    for (int ibody = 0; ibody < static_cast<int>(this->body_size); ++ibody) {
        if (!internally_modified[ibody]) {continue;}
        externally_modified[ibody] = true;
        symmetry_modified[ibody] = std::vector<bool>(this->protein->get_body(ibody).size_symmetry(), true);
    }

    // whether copy isym of a body (0 being the main body) was moved by its symmetry, resp. moved at all
    auto symmetry_changed = [&symmetry_modified] (int ibody, int isym) {return isym != 0 && symmetry_modified[ibody][isym-1];};
    auto moved = [&externally_modified, &symmetry_changed] (int ibody, int isym) {return externally_modified[ibody] || symmetry_changed(ibody, isym);};
    auto copies = [this] (int ibody) {return 1+this->protein->get_body(ibody).size_symmetry();};

    for (int ibody1 = 0; ibody1 < static_cast<int>(this->body_size); ++ibody1) {
        if (internally_modified[ibody1]) {
            calc_aa_self(&calculator, ibody1);
        }

        // a pair of copies of different bodies must be recalculated if either of them moved
        for (int ibody2 = 0; ibody2 < ibody1; ++ibody2) {
            for (int isym1 = 0; isym1 < copies(ibody1); ++isym1) {
                for (int isym2 = 0; isym2 < copies(ibody2); ++isym2) {
                    if (moved(ibody1, isym1) || moved(ibody2, isym2)) {
                        calc_aa(&calculator, ibody1, isym1, ibody2, isym2);
                    }
                }
            }
        }

        // the copies of a body move along with it, so the pairs within it only change with its symmetries
        for (int isym1 = 1; isym1 < copies(ibody1); ++isym1) {
            for (int isym2 = 0; isym2 < isym1; ++isym2) {
                if (symmetry_changed(ibody1, isym1) || symmetry_changed(ibody1, isym2)) {
                    calc_aa(&calculator, ibody1, isym1, ibody1, isym2);
                }
            }
        }

        if constexpr (hydration_enabled) {
            for (int isym = 0; isym < copies(ibody1); ++isym) {
                if (hydration_modified || moved(ibody1, isym)) {
                    calc_aw(&calculator, ibody1, isym);
                }
            }
        }
    }

    // the recalculated partial histograms replace their old contents in the store, which were taken out of the master histogram as they were queued
    calculator.run();
    for (int id : recalculated) {store->visit(id, [this] (const auto& result) {fold_classes(this->master, result, std::plus<>());});}
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
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::update_compact_representation_body(int ibody) {
    coords[ibody] = symmetry::detail::generate_transformed_data<form_factors>(this->protein->get_body(ibody));
    if constexpr (!form_factors) {
        for (auto& c : coords[ibody].atomic) {
            for (auto& sym : c) {
                hist::detail::SimpleExvModel::apply_simple_excluded_volume(sym, this->protein);
            }
        }
    }
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::update_compact_representation_symmetry(int ibody, int isym) {
    assert(ibody < static_cast<int>(coords.size()) && "update_compact_representation_symmetry: ibody out of range");
    assert(isym > 0 && isym < static_cast<int>(coords[ibody].atomic.size()) && "update_compact_representation_symmetry: isym out of range");
    coords[ibody].atomic[isym] = symmetry::detail::generate_transformed_data<form_factors>(this->protein->get_body(ibody), isym-1).data;
    if constexpr (!form_factors) {
        for (auto& sym : coords[ibody].atomic[isym]) {
            hist::detail::SimpleExvModel::apply_simple_excluded_volume(sym, this->protein);
        }
    }
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::update_compact_representation_water() {
    coords_w = hist::detail::factory::construct_from_waters(this->protein);
}

template<bool weighted_bins, bool form_factors>
std::unique_ptr<ICompositeDistanceHistogram> PartialSymmetryManagerMTBase<weighted_bins, form_factors>::calculate_all() {
    auto total = calculate();
    int bins = total->get_weighted_counts().size();

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions.
    // a result that is never calculated is zero, so each contribution is simply the sum over all of its results
    using Distributions = ManagerDistributions<weighted_bins, form_factors>;
    std::vector<int> aa_ids, aw_ids;
    aa.for_each_id([&aa_ids] (int id) {aa_ids.push_back(id);});
    for (const auto& ids : aw) {aw_ids.insert(aw_ids.end(), ids.begin(), ids.end());}

    Distributions d;
    d.p_aa = sum_results<typename Distributions::aa_t>(*store, aa_ids);
    d.p_aw = sum_results<typename Distributions::aw_t>(*store, aw_ids);
    d.p_ww = store->get_1d(ww);
    d.p_tot = this->master; // NOLINT - intentional slicing
    d.resize(bins);
    return make_histogram(std::move(d), exv_method, protein);
}

template<bool weighted_bins, bool form_factors>
int PartialSymmetryManagerMTBase<weighted_bins, form_factors>::prepare_axis() {
    int required = hist::detail::required_partial_bin_count(*this->protein);
    if (!this->master.empty()) {
        if (required <= this->master.axis.bins) {return this->master.axis.bins;}

        logging::log("PartialSymmetryManagerMT::prepare_axis: structure outgrew its axis; rebuilding");
        this->master = hist::detail::MasterHistogram<weighted_bins>();
        this->statemanager->modified_all();
    }
    return hist::detail::grown_partial_bin_count(required);
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::initialize(int bin_count) {
    Axis axis(0, settings::axes::bin_width*bin_count, bin_count);
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram<weighted_bins>(p_base, axis);

    // one result for every calculated pair of symmetries of every body pair, and one for every symmetry of every body
    std::vector<int> sym_counts(this->body_size);
    for (int ibody = 0; ibody < this->body_size; ++ibody) {sym_counts[ibody] = 1 + this->protein->get_body(ibody).size_symmetry();}
    store = std::make_unique<distance_calculator::HistogramStore<weighted_bins>>(axis.bins, form_factors ? form_factor::get_active_count() : 1);
    aa = detail::SymmetryPairIds(sym_counts, [this] () {return form_factors ? store->allocate_3d() : store->allocate_1d();});
    aw.assign(this->body_size, {});
    for (int ibody = 0; ibody < this->body_size; ++ibody) {
        aw[ibody].resize(sym_counts[ibody]);
        for (int& id : aw[ibody]) {id = form_factors ? store->allocate_2d() : store->allocate_1d();}
    }
    ww = store->allocate_1d();

    for (int ibody = 0; ibody < this->body_size; ++ibody) {
        update_compact_representation_body(ibody); //? unnecessary to update whole body; enough to update main body
    }
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::recalculate(int id) {
    // the result is not written until the calculator runs, so its old contents are still there to be taken out
    store->visit(id, [this] (const auto& result) {fold_classes(this->master, result, std::minus<>());});
    recalculated.push_back(id);
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::propagate_reference_symmetry_modifications(
    const std::vector<bool>& externally_modified,
    const std::vector<bool>& internally_modified,
    std::vector<std::vector<bool>>& symmetry_modified
) const {
    for (int ibody = 0; ibody < this->body_size; ++ibody) {
        const auto& body = this->protein->get_body(ibody);
        for (int isym = 0; isym < body.size_symmetry(); ++isym) {
            // a reference symmetry lives on its primary body; the view bodies hold a view that is
            // skipped here, so each group is processed exactly once via its owning instance
            const auto* ref = dynamic_cast<const symmetry::ReferenceSymmetry*>(body.symmetry().get(isym));
            if (ref == nullptr) {continue;}

            // the group's copies are stale if the shared symmetry itself or any participating body
            // (which feeds the combined centre of mass) has changed. The symmetry already knows its
            // members and the slot it occupies on each, so no body search is needed.
            bool affected = false;
            for (std::size_t k = 0; k < ref->bodies.size(); ++k) {
                int b = ref->bodies[k], slot = ref->slots[k];
                affected = affected || symmetry_modified[b][slot] || externally_modified[b] || internally_modified[b];
            }
            if (!affected) {continue;}

            for (std::size_t k = 0; k < ref->bodies.size(); ++k) {symmetry_modified[ref->bodies[k]][ref->slots[k]] = true;}
        }
    }
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::calc_aa_self(calculator_t calculator, int ibody) {
    const auto& body = protein->get_body(ibody);
    // calculate the self correlation within each body and symmetry, equal to (N_sym+1) * (main body self corr)
    int id = aa.id(ibody, 0, ibody, 0);
    recalculate(id);
    calculator->enqueue_calculate_self(coords[ibody].atomic[0][0], id, 1+body.size_symmetry_total());
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::calc_ww(calculator_t calculator) {
    recalculate(ww);
    calculator->enqueue_calculate_self(coords_w, ww);
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::calc_aa(calculator_t calculator, int ibody1, int isym1, int ibody2, int isym2) {
    // every job below accumulates into the same result, so each loop is held and dispatched as a single group.
    int id = aa.id(ibody1, isym1, ibody2, isym2);
    recalculate(id);

    // internal correlations within the same body
    if (ibody1 == ibody2) {
        assert(isym1 != isym2 && "This method is unsuitable for calculating self-correlations");

        // correlations between a symmetry and its host body
        if (isym2 == 0) {
            const auto& body1 = protein->get_body(ibody1);
            assert(isym1 < 1+static_cast<int>(body1.size_symmetry()) && "symmetry index out of bounds");
            const auto* sym1 = body1.symmetry().get(isym1-1);

            // distinct distance pairs among {original, copy_1, ..., copy_N}; repetition 0 is
            // the original body (atomic[0][0]), 1..N are the copies (atomic[isym1][rep-1])
            calculator->hold();
            for (const auto& pair : sym1->internal_pair_schedule()) {
                assert((pair.repA == 0 || pair.repA-1 < static_cast<int>(coords[ibody1].atomic[isym1].size())) && "internal_pair_schedule: repA out of range for atomic copies");
                assert((pair.repB == 0 || pair.repB-1 < static_cast<int>(coords[ibody1].atomic[isym1].size())) && "internal_pair_schedule: repB out of range for atomic copies");
                const auto& atomicA = pair.repA == 0 ? coords[ibody1].atomic[0][0] : coords[ibody1].atomic[isym1][pair.repA-1];
                const auto& atomicB = pair.repB == 0 ? coords[ibody1].atomic[0][0] : coords[ibody1].atomic[isym1][pair.repB-1];
                calculator->enqueue_calculate_cross(atomicA, atomicB, id, 2*pair.scale);
            }
            calculator->release_hold();
            return;
        }
        assert(isym1 != 0 && "Attempting to calculate cross-correlations outside the lower triangle");
    }

    // every copy of the one against every copy of the other; the main body (0) is its own single copy
    assert(isym1 < static_cast<int>(coords[ibody1].atomic.size()) && isym2 < static_cast<int>(coords[ibody2].atomic.size()) && "symmetry index out of bounds");
    calculator->hold();
    for (const auto& copy1 : coords[ibody1].atomic[isym1]) {
        for (const auto& copy2 : coords[ibody2].atomic[isym2]) {
            calculator->enqueue_calculate_cross(copy1, copy2, id, 2);
        }
    }
    calculator->release_hold();
}

template<bool weighted_bins, bool form_factors>
void PartialSymmetryManagerMTBase<weighted_bins, form_factors>::calc_aw(calculator_t calculator, int ibody, int isym) {
    assert(isym < static_cast<int>(aw[ibody].size()) && "PartialSymmetryManagerMT::calc_aw: symmetry index out of range; symmetries may not be added after the first calculation");
    int id = aw[ibody][isym];
    recalculate(id);

    // every copy against the hydration layer; the main body (0) is its own single copy
    calculator->hold(); // one result for every copy, see calc_aa
    for (const auto& copy : coords[ibody].atomic[isym]) {
        calculator->enqueue_calculate_cross(copy, coords_w, id, form_factors ? 1 : 2); // see HistogramManagerMTBase
    }
    calculator->release_hold();
}

template class hist::PartialSymmetryManagerMTBase<false, false>;
template class hist::PartialSymmetryManagerMTBase<false, true>;
template class hist::PartialSymmetryManagerMTBase<true, false>;
template class hist::PartialSymmetryManagerMTBase<true, true>;