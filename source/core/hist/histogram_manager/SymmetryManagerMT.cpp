// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/SymmetryManagerMT.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <form_factor/FormFactorType.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/histogram_manager/detail/ManagerResults.h>
#include <hist/histogram_manager/detail/SymmetryHelpers.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <utility/Logging.h>

#include <ranges>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::hist::detail;
using namespace ausaxs::symmetry::detail;

template<bool weighted_bins, bool form_factors>
hist::SymmetryManagerMTBase<weighted_bins, form_factors>::SymmetryManagerMTBase(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method) 
    : protein(protein), exv_method(exv_method) 
{}

template<bool weighted_bins, bool form_factors>
std::unique_ptr<hist::DistanceHistogram> hist::SymmetryManagerMTBase<weighted_bins, form_factors>::calculate() {
    return calculate_all();
}

template<bool weighted_bins, bool form_factors>
std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMTBase<weighted_bins, form_factors>::calculate_all() {
    if (protein->size_water() == 0) {
        return calculate<false>();
    }
    return calculate<true>();
}

template<bool weighted_bins, bool form_factors> template <bool contains_waters>
std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMTBase<weighted_bins, form_factors>::calculate() {
    logging::log("SymmetryManagerMT::calculate: starting calculation");

    // start by generating the transformed data
    // note that we are responsible for guaranteeing their lifetime until all enqueue_calculate_* calls are done
    auto[data, data_w] = generate_transformed_data<form_factors>(*protein);
    if constexpr (!form_factors) {
        for (auto& body : data) {
            for (auto& copies : body.atomic) {
                for (auto& copy : copies) {hist::detail::SimpleExvModel::apply_simple_excluded_volume(copy, protein);}
            }
        }
    }

    // the per-body data is a struct rather than a range, so project out the coordinate sets for the estimator
    auto atomic = data | std::views::transform([] (const auto& body) -> const auto& {return body.atomic;});
    int bin_count = hist::detail::required_bin_count(atomic, data_w);

    // every self and cross contribution of a kind sums into the same result; with form factors the atoms are partitioned by type
    hist::distance_calculator::HistogramStore<weighted_bins> store(bin_count, form_factors ? form_factor::get_active_count() : 1);
    int aa = form_factors ? store.allocate_3d() : store.allocate_1d();
    int aw = form_factors ? store.allocate_2d() : store.allocate_1d();
    int ww = store.allocate_1d();
    hist::distance_calculator::Calculator<weighted_bins, form_factors> calculator(store);
    constexpr int aw_pair_factor = form_factors ? 1 : 2; // see HistogramManagerMTBase

    const auto& waters = data_w;

    // resolve a (body, symmetry, repetition) triple to its transformed coordinates;
    // repetition 0 is the original body, 1..N are the generated copies
    auto atomic_at = [&data](int i_body, int i_sym, int rep) -> const AtomicCoordinates<form_factors>& {
        return rep == 0 ? data[i_body].atomic[0][0] : data[i_body].atomic[1+i_sym][rep-1];
    };

    for (int i_body1 = 0; i_body1 < protein->size_body(); ++i_body1) {
        const auto& body = protein->get_body(i_body1);
        const auto& body1_atomic = data[i_body1].atomic[0][0];
        // every copy has identical internal distances, so evaluate once and scale
        calculator.enqueue_calculate_self(body1_atomic, aa, 1 + body.size_symmetry_total());
        if constexpr (contains_waters) {
            calculator.enqueue_calculate_cross(body1_atomic, waters, aw, aw_pair_factor);
        }

        for (int i_sym1 = 0; i_sym1 < body.size_symmetry(); ++i_sym1) {
            const auto* sym1 = body.symmetry().get(i_sym1);

            // distinct distance pairs among {original, copy_1, ..., copy_N} of this symmetry;
            // every other copy-pair is identical to a listed representative and folded into scale
            calculator.hold();
            for (const auto& pair : sym1->internal_pair_schedule()) {
                calculator.enqueue_calculate_cross(
                    atomic_at(i_body1, i_sym1, pair.repA),
                    atomic_at(i_body1, i_sym1, pair.repB),
                    aa, 2*pair.scale
                );
            }
            calculator.release_hold();

            for (int i_repeat1 = 0; i_repeat1 < sym1->repetitions(); ++i_repeat1) {
                const auto& body1_sym_atomic = data[i_body1].atomic[1+i_sym1][i_repeat1];

                // this copy against everything it can pair with, as one group
                calculator.hold();
                if constexpr (contains_waters) {
                    calculator.enqueue_calculate_cross(body1_sym_atomic, waters, aw, aw_pair_factor);
                }

                // external histograms with other bodies
                for (int j_body1 = i_body1+1; j_body1 < protein->size_body(); ++j_body1) {
                    const auto& body2 = protein->get_body(j_body1);
                    const auto& body2_atomic = data[j_body1].atomic[0][0];
                    calculator.enqueue_calculate_cross(body2_atomic, body1_sym_atomic, aa, 2);

                    // external histograms with other symmetries in same body
                    for (int j_sym1 = 0; j_sym1 < body2.size_symmetry(); ++j_sym1) {
                        const auto& sym2 = body2.symmetry().get(j_sym1);
                        for (int j_repeat1 = 0; j_repeat1 < sym2->repetitions(); ++j_repeat1) {
                            const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                            calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, aa, 2);
                        }
                    }
                }

                // internal histogram with other symmetries in same body
                for (int i_sym2 = i_sym1+1; i_sym2 < body.size_symmetry(); ++i_sym2) {
                    const auto& sym2 = body.symmetry().get(i_sym2);
                    for (int i_repeat2 = 0; i_repeat2 < sym2->repetitions(); ++i_repeat2) {
                        const auto& body2_sym_atomic = data[i_body1].atomic[1+i_sym2][i_repeat2];
                        calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, aa, 2);
                    }
                }
                calculator.release_hold();
            }
        }

        // external histograms with other bodies
        for (int j_body1 = i_body1+1; j_body1 < protein->size_body(); ++j_body1) {
            const auto& body2 = protein->get_body(j_body1);
            const auto& body2_atomic = data[j_body1].atomic[0][0];

            // the host body against all of body2, as one group
            calculator.hold();
            calculator.enqueue_calculate_cross(body1_atomic, body2_atomic, aa, 2);

            // external histograms with other symmetries in same body
            for (int j_sym1 = 0; j_sym1 < body2.size_symmetry(); ++j_sym1) {
                const auto& sym2 = body2.symmetry().get(j_sym1);
                for (int j_repeat1 = 0; j_repeat1 < sym2->repetitions(); ++j_repeat1) {
                    const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                    calculator.enqueue_calculate_cross(body1_atomic, body2_sym_atomic, aa, 2);
                }
            }
            calculator.release_hold();
        }
    }
    if constexpr (contains_waters) {
        calculator.enqueue_calculate_self(waters, ww);
    }
    calculator.run();

    // without waters, aw and ww were never named, and are still zero
    return hist::detail::make_histogram(hist::detail::export_distributions<weighted_bins, form_factors>(store, aa, aw, ww), exv_method, protein);
}

template class hist::SymmetryManagerMTBase<false, false>;
template class hist::SymmetryManagerMTBase<false, true>;
template class hist::SymmetryManagerMTBase<true, false>;
template class hist::SymmetryManagerMTBase<true, true>;
