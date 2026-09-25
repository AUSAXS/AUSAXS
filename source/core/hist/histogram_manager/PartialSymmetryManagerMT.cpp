// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/PartialSymmetryManagerMT.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/state/StateManager.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/histogram_manager/detail/PartialBinEstimate.h>
#include <hist/histogram_manager/detail/SymmetryHelpers.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
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

template<bool weighted_bins, bool variable_bin_width> 
PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::PartialSymmetryManagerMT(observer_ptr<const data::Molecule> protein) 
    : IPartialHistogramManager(protein), 
      protein(protein),
      coords(this->body_size)
{}

template<bool weighted_bins, bool variable_bin_width> 
PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::~PartialSymmetryManagerMT() = default;

template<bool weighted_bins, bool variable_bin_width>
std::unique_ptr<DistanceHistogram> PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::calculate() {
    logging::log("PartialSymmetryManagerMT::calculate: starting calculation");
    if (protein->size_water() == 0 && !this->statemanager->is_modified_hydration()) {
        return _calculate<false>();
    }
    return _calculate<true>();
}

template<bool weighted_bins, bool variable_bin_width> template<bool hydration_enabled>
std::unique_ptr<DistanceHistogram> PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::_calculate() {
    if (!this->statemanager->is_modified() && !cache.p_tot.empty()) {
        logging::log("PartialSymmetryManagerMT::calculate: returning cached value");
        auto p_tot = cache.p_tot; // if the state was not modified, we can return the cached value
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

    distance_calculator::Calculator<weighted_bins, variable_bin_width> calculator(*store);

    if constexpr (hydration_enabled) {
        // check if the hydration layer was modified
        if (hydration_modified) {
            calc_ww(&calculator);
        }
    }

    // iterate through the lower triangle and check if either of each pair of different bodies were modified
    for (int ibody1 = 0; ibody1 < static_cast<int>(this->body_size); ++ibody1) {
        // check for internal modifications
        if (internally_modified[ibody1]) {
            calc_aa_self(&calculator, ibody1);

            // internal modification implies external modification
            // everything connected to this body must be recalculated
            externally_modified[ibody1] = true;
            symmetry_modified[ibody1] = std::vector<bool>(this->protein->get_body(ibody1).size_symmetry(), true);
        }

        // check for external modifications 
        for (int ibody2 = 0; ibody2 < ibody1; ++ibody2) {
            // note: off-diagonal elements only
            if (externally_modified[ibody1] || externally_modified[ibody2]) {

                // external modification requires recalculation of all affected symmetries
                for (int isym1 = 0; isym1 < 1+this->protein->get_body(ibody1).size_symmetry(); ++isym1) {
                    for (int isym2 = 0; isym2 < 1+this->protein->get_body(ibody2).size_symmetry(); ++isym2) {
                        calc_aa(&calculator, ibody1, isym1, ibody2, isym2);
                    }
                }
            }
        }

        // if everything was not already recalculated due to an external modification in the previous loop, 
        // we must check for modifications to each symmetry
        if (!externally_modified[ibody1]) {

            // correlations between this main body and symmetries of other main bodies
            for (int ibody2 = 0; ibody2 < ibody1; ++ibody2) {
                if (externally_modified[ibody2]) {continue;}
                for (int isym2 = 0; isym2 < this->protein->get_body(ibody2).size_symmetry(); ++isym2) {
                    if (symmetry_modified[ibody2][isym2]) {
                        calc_aa(&calculator, ibody1, 0, ibody2, isym2+1);
                    }
                }
            }

            // correlations between symmetries of this body and other bodies
            for (int ibody2 = 0; ibody2 < ibody1; ++ibody2) {
                if (externally_modified[ibody2]) {continue;} // already handled
                for (int isym1 = 0; isym1 < this->protein->get_body(ibody1).size_symmetry(); ++isym1) {

                    // cross-correlation with other main body
                    if (symmetry_modified[ibody1][isym1]) {
                        calc_aa(&calculator, ibody1, isym1+1, ibody2, 0);
                    }
                    
                    // cross-correlations with symmetries in other main body
                    for (int isym2 = 0; isym2 < this->protein->get_body(ibody2).size_symmetry(); ++isym2) {
                        if (!(symmetry_modified[ibody1][isym1] || symmetry_modified[ibody2][isym2])) {continue;}
                        calc_aa(&calculator, ibody1, isym1+1, ibody2, isym2+1);
                    }
                }
            }
        }

        // diagonal elements only have to be recalculated if the symmetry was modified
        {
            for (int isym1 = 0; isym1 < this->protein->get_body(ibody1).size_symmetry(); ++isym1) {
                // cross-correlation with main body
                if (symmetry_modified[ibody1][isym1]) {
                    calc_aa(&calculator, ibody1, isym1+1, ibody1, 0);
                }

                // cross-correlations with other symmetries
                for (int isym2 = 0; isym2 < isym1; ++isym2) {
                    if (!(symmetry_modified[ibody1][isym1] || symmetry_modified[ibody1][isym2])) {continue;}
                    calc_aa(&calculator, ibody1, isym1+1, ibody1, isym2+1);
                }
            }
        }

        if constexpr (hydration_enabled) {
            // we also have to remember to update the partial histograms with the hydration layer
            if (externally_modified[ibody1] || hydration_modified) {

                // update all by looping from 0 (main body) to 1+size_symmetry (last symmetry)
                for (int isym1 = 0; isym1 < 1+this->protein->get_body(ibody1).size_symmetry(); ++isym1) {
                    calc_aw(&calculator, ibody1, isym1);
                }
            } else {

                // hydration layer not modified, check for symmetry modifications
                for (int isym1 = 0; isym1 < this->protein->get_body(ibody1).size_symmetry(); ++isym1) {
                    if (symmetry_modified[ibody1][isym1]) {
                        calc_aw(&calculator, ibody1, isym1+1);
                    }
                }
            }
        }
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
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::update_compact_representation_body(int ibody) {
    coords[ibody] = symmetry::detail::generate_transformed_data<variable_bin_width>(this->protein->get_body(ibody));
    for (auto& c : coords[ibody].atomic) {
        for (auto& sym : c) {
            hist::detail::SimpleExvModel::apply_simple_excluded_volume(sym, this->protein);
        }
    }
}

template<bool weighted_bins, bool variable_bin_width>
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::update_compact_representation_symmetry(int ibody, int isym) {
    assert(ibody < static_cast<int>(coords.size()) && "update_compact_representation_symmetry: ibody out of range");
    assert(isym > 0 && isym < static_cast<int>(coords[ibody].atomic.size()) && "update_compact_representation_symmetry: isym out of range");
    coords[ibody].atomic[isym] = symmetry::detail::generate_transformed_data<variable_bin_width>(this->protein->get_body(ibody), isym-1).data;
    for (auto& sym : coords[ibody].atomic[isym]) {
        hist::detail::SimpleExvModel::apply_simple_excluded_volume(sym, this->protein);
    }
}

template<bool weighted_bins, bool variable_bin_width>
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::update_compact_representation_water() {
    coords_w = hist::detail::factory::construct_from_waters<variable_bin_width>(this->protein);
}

template<bool weighted_bins, bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::calculate_all() {
    if (
        !this->statemanager->is_modified() 
        && !cache.p_tot.empty() && !cache.p_aa.empty() && !cache.p_aw.empty() && !cache.p_ww.empty()
    ) {
        logging::log("PartialSymmetryManagerMT::calculate_all: returning cached value");
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

    auto total = calculate();
    int bins = total->get_weighted_counts().size();

    // determine p_tot
    GenericDistribution1D_t p_tot(bins);
    for (int i = 0; i < bins; ++i) {
        p_tot.index(i) = this->master.index(i);
    }

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions.
    // the aa rows come first in the store, then the aw rows, then ww; a row that is never calculated is zero, so each
    // contribution is simply the sum over its block of rows
    GenericDistribution1D_t p_ww = store->export_1d(handle_ww());
    GenericDistribution1D_t p_aa = this->master.base;
    GenericDistribution1D_t p_aw(bins);
    p_ww.resize(bins);
    p_aa.resize(bins);
    for (int h = 0; h < handle_aw(0, 0); ++h) {
        std::transform(p_aa.begin(), p_aa.end(), store->row(h).begin(), p_aa.begin(), std::plus<>());
    }
    for (int h = handle_aw(0, 0); h < handle_ww(); ++h) {
        std::transform(p_aw.begin(), p_aw.end(), store->row(h).begin(), p_aw.begin(), std::plus<>());
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
int PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::prepare_axis() {
    int required = hist::detail::required_partial_bin_count<variable_bin_width>(*this->protein);
    if (!this->master.empty()) {
        if (required <= this->master.axis.bins) {return this->master.axis.bins;}

        logging::log("PartialSymmetryManagerMT::prepare_axis: structure outgrew its axis; rebuilding");
        this->master = hist::detail::MasterHistogram<weighted_bins>();
        this->statemanager->modified_all();
    }
    return hist::detail::grown_partial_bin_count(required);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::initialize(int bin_count) {
    Axis axis(0, settings::axes::bin_width*bin_count, bin_count);
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram<weighted_bins>(p_base, axis);

    // lay out the rows: the aa partials of every body pair, then the aw partials of every body, then ww. a pair of different
    // bodies holds every combination of their symmetries, while a body with itself only holds a lower triangle, as the
    // upper one would duplicate it
    sym_count.resize(this->body_size);
    aa_offset.resize(this->body_size*(this->body_size+1)/2);
    aw_offset.resize(this->body_size);
    int rows = 0;
    for (int ibody = 0; ibody < this->body_size; ++ibody) {
        sym_count[ibody] = 1 + this->protein->get_body(ibody).size_symmetry();
        for (int ibody2 = 0; ibody2 <= ibody; ++ibody2) {
            aa_offset[ibody*(ibody+1)/2 + ibody2] = rows;
            rows += ibody2 == ibody ? sym_count[ibody]*(sym_count[ibody]+1)/2 : sym_count[ibody]*sym_count[ibody2];
        }
    }
    for (int ibody = 0; ibody < this->body_size; ++ibody) {
        aw_offset[ibody] = rows;
        rows += sym_count[ibody];
    }
    store = std::make_unique<distance_calculator::HistogramStore<weighted_bins>>(rows+1, axis.bins); // +1 for ww

    for (int ibody = 0; ibody < this->body_size; ++ibody) {
        update_compact_representation_body(ibody); //? unnecessary to update whole body; enough to update main body
    }
}

template<bool weighted_bins, bool variable_bin_width>
int PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::handle_aa(int ibody1, int isym1, int ibody2, int isym2) const {
    assert(0 <= ibody2 && ibody2 <= ibody1 && ibody1 < this->body_size && "PartialSymmetryManagerMT::handle_aa: expected a body pair in the lower triangle");
    assert(0 <= isym1 && isym1 < sym_count[ibody1] && 0 <= isym2 && isym2 < sym_count[ibody2] && "PartialSymmetryManagerMT::handle_aa: symmetry index out of range; symmetries may not be added after the first calculation");
    int offset = aa_offset[ibody1*(ibody1+1)/2 + ibody2];
    if (ibody1 == ibody2) {
        assert(isym2 <= isym1 && "PartialSymmetryManagerMT::handle_aa: expected a symmetry pair in the lower triangle");
        return offset + isym1*(isym1+1)/2 + isym2;
    }
    return offset + isym1*sym_count[ibody2] + isym2;
}

template<bool weighted_bins, bool variable_bin_width>
int PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::handle_aw(int ibody, int isym) const {
    assert(0 <= isym && isym < sym_count[ibody] && "PartialSymmetryManagerMT::handle_aw: symmetry index out of range; symmetries may not be added after the first calculation");
    return aw_offset[ibody] + isym;
}

template<bool weighted_bins, bool variable_bin_width>
int PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::handle_ww() const {
    return store->rows()-1;
}

template<bool weighted_bins, bool variable_bin_width>
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::recalculate(int h) {
    // the row is not written until the calculator runs, so its old contents are still there to be taken out
    this->master -= store->row(h);
    recalculated.push_back(h);
}

template<bool weighted_bins, bool variable_bin_width>
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::propagate_reference_symmetry_modifications(
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

template<bool weighted_bins, bool variable_bin_width>
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::calc_aa_self(calculator_t calculator, int ibody) {
    const auto& body = protein->get_body(ibody);
    // calculate the self correlation within each body and symmetry, equal to (N_sym+1) * (main body self corr)
    int h = handle_aa(ibody, 0, ibody, 0);
    recalculate(h);
    calculator->enqueue_calculate_self(coords[ibody].atomic[0][0], h, 1+body.size_symmetry_total());
}

template<bool weighted_bins, bool variable_bin_width> 
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::calc_ww(calculator_t calculator) {
    recalculate(handle_ww());
    calculator->enqueue_calculate_self(coords_w, handle_ww());
}

template<bool weighted_bins, bool variable_bin_width> 
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::calc_aa(calculator_t calculator, int ibody1, int isym1, int ibody2, int isym2) {
    // every job below accumulates into the same row, so each loop is held and dispatched as a single group.
    const auto& body1 = protein->get_body(ibody1);
    const auto& body2 = protein->get_body(ibody2);
    int h = handle_aa(ibody1, isym1, ibody2, isym2);
    recalculate(h);

    // internal correlations within the same body
    if (ibody1 == ibody2) {
        assert(isym1 != isym2 && "This method is unsuitable for calculating self-correlations");

        // correlations between a symmetry and its host body
        if (isym2 == 0) {
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
                calculator->enqueue_calculate_cross(atomicA, atomicB, h, 2*pair.scale);
            }
            calculator->release_hold();
            return;
        }
        assert(isym1 != 0 && "Attempting to calculate cross-correlations outside the lower triangle");
    }

    // symmetry 0 is the main body, so we have to treat it separately
    if (isym1 == 0 && isym2 == 0) {
        calculator->enqueue_calculate_cross(coords[ibody1].atomic[0][0], coords[ibody2].atomic[0][0], h, 2);
        return;
    }
    if (isym1 == 0) {
        assert(isym2 < 1+body2.size_symmetry() && "symmetry index out of bounds");
        const auto& sym2 = body2.symmetry().get(isym2-1);
        calculator->hold();
        for (int irepeat2 = 0; irepeat2 < sym2->repetitions(); ++irepeat2) {
            const auto& body2_sym_atomic = coords[ibody2].atomic[isym2][irepeat2];
            calculator->enqueue_calculate_cross(coords[ibody1].atomic[0][0], body2_sym_atomic, h, 2);
        }
        calculator->release_hold();
        return;
    }
    if (isym2 == 0) {
        assert(isym1 < 1+body1.size_symmetry() && "symmetry index out of bounds");
        const auto& sym1 = body1.symmetry().get(isym1-1);
        calculator->hold();
        for (int irepeat1 = 0; irepeat1 < sym1->repetitions(); ++irepeat1) {
            const auto& body1_sym_atomic = coords[ibody1].atomic[isym1][irepeat1];
            calculator->enqueue_calculate_cross(body1_sym_atomic, coords[ibody2].atomic[0][0], h, 2);
        }
        calculator->release_hold();
        return;
    }

    // else iterate over the replications of both symmetries
    assert(isym1 < 1+body1.size_symmetry() && "symmetry index out of bounds");
    assert(isym2 < 1+body2.size_symmetry() && "symmetry index out of bounds");
    const auto& sym1 = body1.symmetry().get(isym1-1);
    const auto& sym2 = body2.symmetry().get(isym2-1);

    calculator->hold();
    for (int irepeat1 = 0; irepeat1 < sym1->repetitions(); ++irepeat1) {
        const auto& body1_sym_atomic = coords[ibody1].atomic[isym1][irepeat1];
        for (int irepeat2 = 0; irepeat2 < sym2->repetitions(); ++irepeat2) {
            const auto& body2_sym_atomic = coords[ibody2].atomic[isym2][irepeat2];
            calculator->enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, h, 2);
        }
    }
    calculator->release_hold();
}

template<bool weighted_bins, bool variable_bin_width> 
void PartialSymmetryManagerMT<weighted_bins, variable_bin_width>::calc_aw(calculator_t calculator, int ibody, int isym) {
    const auto& body = protein->get_body(ibody);
    const auto& waters = coords_w;
    int h = handle_aw(ibody, isym);
    recalculate(h);

    // symmetry 0 is the main body, so we have to treat it separately
    if (isym == 0) {
        calculator->enqueue_calculate_cross(coords[ibody].atomic[0][0], waters, h, 2);
        return;
    }

    // else iterate over its repititions
    assert(isym < 1+static_cast<int>(body.size_symmetry()) && "symmetry index out of bounds");
    const auto& sym = body.symmetry().get(isym-1);
    calculator->hold(); // one row for every copy, see calc_aa
    for (int irepeat = 0; irepeat < sym->repetitions(); ++irepeat) {
        const auto& body1_sym_atomic = coords[ibody].atomic[isym][irepeat];
        calculator->enqueue_calculate_cross(body1_sym_atomic, waters, h, 2);
    }
    calculator->release_hold();
}

template class hist::PartialSymmetryManagerMT<false, false>;
template class hist::PartialSymmetryManagerMT<false, true>;
template class hist::PartialSymmetryManagerMT<true, false>;
template class hist::PartialSymmetryManagerMT<true, true>;