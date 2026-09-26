// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/SymmetryManagerMT.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/histogram_manager/detail/SymmetryHelpers.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <utility/Logging.h>

#include <ranges>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::hist::detail;
using namespace ausaxs::symmetry::detail;

template<bool weighted_bins>
hist::SymmetryManagerMT<weighted_bins>::SymmetryManagerMT(observer_ptr<const data::Molecule> protein) : protein(protein) {}

template<bool weighted_bins>
std::unique_ptr<hist::DistanceHistogram> hist::SymmetryManagerMT<weighted_bins>::calculate() {
    return calculate_all();
}

template<bool weighted_bins>
std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<weighted_bins>::calculate_all() {
    if (protein->size_water() == 0) {
        return calculate<false>();
    }
    return calculate<true>();
}

template<bool weighted_bins> template <bool contains_waters>
std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<weighted_bins>::calculate() {
    logging::log("SymmetryManagerMT::calculate: starting calculation");

    using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;

    // start by generating the transformed data
    // note that we are responsible for guaranteeing their lifetime until all enqueue_calculate_* calls are done
    auto[data, data_w] = generate_transformed_data(*protein);

    // the per-body data is a struct rather than a range, so project out the coordinate sets for the estimator
    auto atomic = data | std::views::transform([] (const auto& body) -> const auto& {return body.atomic;});
    int bin_count = hist::detail::required_bin_count(atomic, data_w);

    // every self and cross contribution of a kind sums into the same row
    hist::distance_calculator::HistogramStore<weighted_bins> store(bin_count);
    int aa = store.allocate_1d(), aw = store.allocate_1d(), ww = store.allocate_1d();
    hist::distance_calculator::Calculator<weighted_bins> calculator(store);

    const auto& waters = data_w;

    // resolve a (body, symmetry, repetition) triple to its transformed coordinates;
    // repetition 0 is the original body, 1..N are the generated copies
    auto atomic_at = [&data](int i_body, int i_sym, int rep) -> const CompactCoordinates& {
        return rep == 0 ? data[i_body].atomic[0][0] : data[i_body].atomic[1+i_sym][rep-1];
    };

    for (int i_body1 = 0; i_body1 < protein->size_body(); ++i_body1) {
        const auto& body = protein->get_body(i_body1);
        const auto& body1_atomic = data[i_body1].atomic[0][0];
        // every copy has identical internal distances, so evaluate once and scale
        calculator.enqueue_calculate_self(body1_atomic, aa, 1 + body.size_symmetry_total());
        if constexpr (contains_waters) {
            calculator.enqueue_calculate_cross(waters, body1_atomic, aw, 2);
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
                    calculator.enqueue_calculate_cross(waters, body1_sym_atomic, aw, 2);
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
    GenericDistribution1D_t p_aa = store.export_1d(aa);
    GenericDistribution1D_t p_aw = store.export_1d(aw);
    GenericDistribution1D_t p_ww = store.export_1d(ww);

    // calculate p_tot
    GenericDistribution1D_t p_tot(bin_count);
    for (int i = 0; i < static_cast<int>(p_tot.size()); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + p_aw.index(i);}

    // downsize our axes to only the relevant area
    int max_bin = hist::detail::trimmed_bin_count(p_tot);
    p_aa.resize(max_bin);
    p_ww.resize(max_bin);
    p_aw.resize(max_bin);
    p_tot.resize(max_bin);

    if constexpr (weighted_bins) {
        return std::make_unique<hist::CompositeDistanceHistogram>(
            hist::Distribution1D(std::move(p_aa)), 
            hist::Distribution1D(std::move(p_aw)), 
            hist::Distribution1D(std::move(p_ww)), 
            std::move(p_tot)
        );
    } else {
        return std::make_unique<hist::CompositeDistanceHistogram>(
            std::move(p_aa), 
            std::move(p_aw), 
            std::move(p_ww), 
            std::move(p_tot)
        );
    }
}

template class hist::SymmetryManagerMT<false>;
template class hist::SymmetryManagerMT<true>;
