#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/detail/SymmetryHelpers.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/detail/SimpleExvModel.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::hist::detail;
using namespace ausaxs::symmetry::detail;

template<bool use_weighted_distribution>
hist::SymmetryManagerMT<use_weighted_distribution>::SymmetryManagerMT(observer_ptr<const data::Molecule> protein) : protein(protein) {}

template<bool use_weighted_distribution>
std::unique_ptr<hist::DistanceHistogram> hist::SymmetryManagerMT<use_weighted_distribution>::calculate() {
    return calculate_all();
}

template<bool use_weighted_distribution>
std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<use_weighted_distribution>::calculate_all() {
    if (protein->size_water() == 0) {
        return calculate<false>();
    } else {
        return calculate<true>();
    }
}

template<bool use_weighted_distribution> template <bool contains_waters>
std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<use_weighted_distribution>::calculate() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    hist::distance_calculator::SimpleCalculator<use_weighted_distribution> calculator;

    // start by generating the transformed data
    // note that we are responsible for guaranteeing their lifetime until all enqueue_calculate_* calls are done
    auto[data, data_w] = generate_transformed_data(*protein);

    const auto& waters = data_w;
    int self_merge_id_aa = 0, self_merge_id_ww = 1;
    int cross_merge_id_aa = 0, cross_merge_id_aw = 1, cross_merge_id_ww = 2;
    for (int i_body1 = 0; i_body1 < static_cast<int>(protein->size_body()); ++i_body1) {
        const auto& body = protein->get_body(i_body1);
        const auto& body1_atomic = data[i_body1].atomic[0][0];
        calculator.enqueue_calculate_self(body1_atomic, 1 + body.size_symmetry_total(), self_merge_id_aa);
        if constexpr (contains_waters) {
            calculator.enqueue_calculate_cross(waters, body1_atomic, 1, cross_merge_id_aw);
        }

        for (int i_sym1 = 0; i_sym1 < static_cast<int>(body.size_symmetry()); ++i_sym1) {
            const auto& sym1 = body.symmetry().get(i_sym1);
            bool closed = sym1.is_closed() && sym1.repeat_relation.rotate.magnitude() != 0;
            for (int i_repeat1 = 0; i_repeat1 < (sym1.repetitions - closed); ++i_repeat1) {
                const auto& body1_sym_atomic = data[i_body1].atomic[1+i_sym1][i_repeat1];

                // assume we have 3 repeats of symmetry B, so we have the bodies: A B1 B2 B3. Then
                // AB1 == AB2 == AB3, (scale AB1 by 3)
                // AB2 == B1B3        (scale B1B3 by 2)
                // AB3                (scale AB3 by 1)
                //
                // if the symmetry is closed, AB3 == AB1
                int scale = sym1.repetitions - i_repeat1;
                if (i_repeat1 == 0 && closed) {scale += 1;}
                calculator.enqueue_calculate_cross(body1_atomic, body1_sym_atomic, scale, cross_merge_id_aa);
            }

            for (int i_repeat1 = 0; i_repeat1 < sym1.repetitions; ++i_repeat1) {
                const auto& body1_sym_atomic = data[i_body1].atomic[1+i_sym1][i_repeat1];
                if constexpr (contains_waters) {
                    calculator.enqueue_calculate_cross(waters, body1_sym_atomic, 1, cross_merge_id_aw);
                }

                // external histograms with other bodies
                for (int j_body1 = i_body1+1; j_body1 < static_cast<int>(protein->size_body()); ++j_body1) {
                    const auto& body2 = protein->get_body(j_body1);
                    const auto& body2_atomic = data[j_body1].atomic[0][0];
                    calculator.enqueue_calculate_cross(body2_atomic, body1_sym_atomic, 1, cross_merge_id_aa);

                    // external histograms with other symmetries in same body
                    for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                        const auto& sym2 = body2.symmetry().get(j_sym1);
                        for (int j_repeat1 = 0; j_repeat1 < sym2.repetitions; ++j_repeat1) {
                            const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                            calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, 1, cross_merge_id_aa);
                        }
                    }
                }

                // internal histogram with other symmetries in same body
                for (int i_sym2 = i_sym1+1; i_sym2 < static_cast<int>(body.size_symmetry()); ++i_sym2) {
                    const auto& sym2 = body.symmetry().get(i_sym2);
                    for (int i_repeat2 = 0; i_repeat2 < sym2.repetitions; ++i_repeat2) {
                        const auto& body2_sym_atomic = data[i_body1].atomic[1+i_sym2][i_repeat2];
                        calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, 1, cross_merge_id_aa);
                    }
                }
            }
        }

        // external histograms with other bodies
        for (int j_body1 = i_body1+1; j_body1 < static_cast<int>(protein->size_body()); ++j_body1) {
            const auto& body2 = protein->get_body(j_body1);
            const auto& body2_atomic = data[j_body1].atomic[0][0];
            calculator.enqueue_calculate_cross(body1_atomic, body2_atomic, 1, cross_merge_id_aa);

            // external histograms with other symmetries in same body
            for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                const auto& sym2 = body2.symmetry().get(j_sym1);
                for (int j_repeat1 = 0; j_repeat1 < sym2.repetitions; ++j_repeat1) {
                    const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                    calculator.enqueue_calculate_cross(body1_atomic, body2_sym_atomic, 1, cross_merge_id_aa);
                }
            }
        }
    }
    if constexpr (contains_waters) {
        calculator.enqueue_calculate_self(waters, 1, self_merge_id_ww);
    }

    auto res = calculator.run();
    assert((contains_waters ? 2 : 1) == res.self.size() && "SymmetryManager::calculate: self size mismatch");

    GenericDistribution1D_t p_aa = std::move(res.self[self_merge_id_aa]);
    GenericDistribution1D_t p_ww, p_aw;

    // merge results
    if (res.cross.contains(cross_merge_id_aa)) {
        p_aa += res.cross[ cross_merge_id_aa];
    }

    if constexpr (contains_waters) {
        if (res.cross.contains(        cross_merge_id_aw)) {
            p_aw = std::move(res.cross[cross_merge_id_aw]);
        }
        if (res.self.contains(        self_merge_id_ww)) {
            p_ww = std::move(res.self[self_merge_id_ww]);
        }
        if (res.cross.contains(cross_merge_id_ww)) {
            p_ww += res.cross[ cross_merge_id_ww];
        }
    } else {
        p_ww = GenericDistribution1D_t(p_aa.size());
        p_aw = GenericDistribution1D_t(p_aa.size());
    }

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

template std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<true>::calculate<true>();
template std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<true>::calculate<false>();
template std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<false>::calculate<true>();
template std::unique_ptr<hist::ICompositeDistanceHistogram> hist::SymmetryManagerMT<false>::calculate<false>();
template class hist::SymmetryManagerMT<true>;
template class hist::SymmetryManagerMT<false>;