#include <data/symmetry/SymmetryManagerMT.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/detail/SimpleExvModel.h>

#include <cassert>

using namespace ausaxs;
using namespace hist::detail;

namespace local {
    enum h_type {
        SELF_AA, SELF_AW, SELF_WW, 
        SELF_SYM_AA, SELF_SYM_AW, SELF_SYM_WW,
        CROSS_AA, CROSS_AW, CROSS_WW
    };

    struct BodySymmetryData {
        template<typename T>
        using symmetry_t = std::vector<T>;

        template<typename T>
        using repetition_t = std::vector<T>;

        // the outer loop is over the symmetries, the inner loop is over the repetitions
        // index [0][0] is the original data
        symmetry_t<repetition_t<CompactCoordinates>> atomic;
        CompactCoordinates waters;
    };

    struct ScaleResult {
        ScaleResult() = default;
        ScaleResult(int index, int scale) : index(index), scale(scale) {}
        ScaleResult(int index, std::size_t scale) : index(index), scale(static_cast<int>(scale)) {}
        int index;
        int scale;
    };

    std::vector<BodySymmetryData> generate_transformed_data(const data::Molecule& protein) {
        std::vector<BodySymmetryData> res(protein.size_body());

        // for every body
        for (int i_body1 = 0; i_body1 < static_cast<int>(protein.size_body()); ++i_body1) {
            const auto& body = protein.get_body(i_body1);
            CompactCoordinates data_a(body.get_atoms());
            CompactCoordinates data_w;
            if (0 < protein.size_water()) {
                data_w = CompactCoordinates(body.get_waters());
            }
            hist::detail::SimpleExvModel::apply_simple_excluded_volume(data_a, &protein);

            // loop over its symmetries
            std::vector<std::vector<CompactCoordinates>> atomic(1+body.size_symmetry());
            for (int i_sym_1 = 0; i_sym_1 < static_cast<int>(body.size_symmetry()); ++i_sym_1) {
                const auto& symmetry = body.symmetry().get(i_sym_1);
                std::vector<CompactCoordinates> sym_atomic(symmetry.repeat, data_a);

                // for every symmetry, loop over how many times it should be repeated
                // it is then repeatedly applied to the same data
                for (int i_repeat = 0; i_repeat < symmetry.repeat; ++i_repeat) {
                    auto t = symmetry.get_transform<double>(i_repeat+1);
                    std::transform(
                        sym_atomic[i_repeat].get_data().begin(), 
                        sym_atomic[i_repeat].get_data().end(), 
                        sym_atomic[i_repeat].get_data().begin(), 
                        [t] (const CompactCoordinatesData& v) -> CompactCoordinatesData {return {t(v.value.pos), v.value.w}; }
                    );
                }
                atomic[1+i_sym_1] = std::move(sym_atomic);
            }
            atomic[0]    = {std::move(data_a)};
            res[i_body1] = {std::move(atomic), std::move(data_w)};
        }
        return res;
    }
}

template<bool use_weighted_distribution>
std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate(const data::Molecule& protein) {
    if (protein.size_water() == 0) {
        return calculate<use_weighted_distribution, false>(protein);
    } else {
        return calculate<use_weighted_distribution, true>(protein);
    }
}

template<bool use_weighted_distribution, bool contains_waters>
std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate(const data::Molecule& protein) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    hist::distance_calculator::SimpleCalculator<use_weighted_distribution> calculator;

    // start by generating the transformed data
    // note that we are responsible for guaranteeing their lifetime until all enqueue_calculate_* calls are done
    auto data = local::generate_transformed_data(protein);

    int self_merge_id_aa = 0, self_merge_id_ww = 1;
    int cross_merge_id_aa = 0, cross_merge_id_aw = 1, cross_merge_id_ww = 2;
    for (int i_body1 = 0; i_body1 < static_cast<int>(protein.size_body()); ++i_body1) {
        const auto& body = protein.get_body(i_body1);
        const auto& body1_atomic = data[i_body1].atomic[0][0];
        const auto& body1_waters = data[i_body1].waters;
        calculator.enqueue_calculate_self(body1_atomic, 1 + body.size_symmetry_total(), self_merge_id_aa);
        if constexpr (contains_waters) {
            calculator.enqueue_calculate_self(body1_waters, 1, self_merge_id_ww);
            calculator.enqueue_calculate_cross(body1_atomic, body1_waters, 1, cross_merge_id_aw);
        }

        for (int i_sym1 = 0; i_sym1 < static_cast<int>(body.size_symmetry()); ++i_sym1) {
            const auto& sym1 = body.symmetry().get(i_sym1);
            bool closed = sym1.is_closed();
            for (int i_repeat1 = 0; i_repeat1 < (sym1.repeat - closed); ++i_repeat1) {
                const auto& body1_sym_atomic = data[i_body1].atomic[1+i_sym1][i_repeat1];

                // assume we have 3 repeats of symmetry B, so we have the bodies: A B1 B2 B3. Then
                // AB1 == AB2 == AB3, (scale AB1 by 3)
                // AB2 == B1B3        (scale B1B3 by 2)
                // AB3                (scale AB3 by 1)
                //
                // if the symmetry is closed, AB3 == AB1
                int scale = sym1.repeat - i_repeat1;
                if (i_repeat1 == 0 && closed) {scale += 1;}

                calculator.enqueue_calculate_cross(body1_atomic, body1_sym_atomic, scale, cross_merge_id_aa);
                if constexpr (contains_waters) {
                    calculator.enqueue_calculate_cross(body1_waters, body1_sym_atomic, scale, cross_merge_id_aw);
                }

                // external histograms with other bodies
                for (int j_body1 = i_body1+1; j_body1 < static_cast<int>(protein.size_body()); ++j_body1) {
                    const auto& body2 = protein.get_body(j_body1);
                    const auto& body2_atomic = data[j_body1].atomic[0][0];
                    const auto& body2_waters = data[j_body1].waters;
                    calculator.enqueue_calculate_cross(body2_atomic, body1_sym_atomic, 1, cross_merge_id_aa);
                    if constexpr (contains_waters) {
                        calculator.enqueue_calculate_cross(body2_waters, body1_sym_atomic, 1, cross_merge_id_aw);
                    }

                    // external histograms with other symmetries in same body
                    for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                        const auto& sym2 = body2.symmetry().get(j_sym1);
                        for (int j_repeat1 = 0; j_repeat1 < sym2.repeat; ++j_repeat1) {
                            const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                            calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, 1, cross_merge_id_aa);
                        }
                    }
                }

                // internal histogram with other symmetries in same body
                for (int i_sym2 = i_sym1+1; i_sym2 < static_cast<int>(body.size_symmetry()); ++i_sym2) {
                    const auto& sym2 = body.symmetry().get(i_sym2);
                    for (int i_repeat2 = 0; i_repeat2 < sym2.repeat; ++i_repeat2) {
                        const auto& body2_sym_atomic = data[i_body1].atomic[1+i_sym2][i_repeat2];
                        calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, 1, cross_merge_id_aa);
                    }
                }
            }
        }

        // external histograms with other bodies
        for (int j_body1 = i_body1+1; j_body1 < static_cast<int>(protein.size_body()); ++j_body1) {
            const auto& body2 = protein.get_body(j_body1);
            const auto& body2_atomic = data[j_body1].atomic[0][0];
            const auto& body2_waters = data[j_body1].waters;

            calculator.enqueue_calculate_cross(body1_atomic, body2_atomic, 1, cross_merge_id_aa);
            if constexpr (contains_waters) {
                calculator.enqueue_calculate_cross(body1_atomic, body2_waters, 1, cross_merge_id_aw);
                calculator.enqueue_calculate_cross(body1_waters, body2_waters, 1, cross_merge_id_ww);
                calculator.enqueue_calculate_cross(body1_waters, body2_atomic, 1, cross_merge_id_aw);
            }

            // external histograms with other symmetries in same body
            for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                const auto& sym2 = body2.symmetry().get(j_sym1);
                for (int j_repeat1 = 0; j_repeat1 < sym2.repeat; ++j_repeat1) {
                    const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];

                    calculator.enqueue_calculate_cross(body1_atomic, body2_sym_atomic, 1, cross_merge_id_aa);
                    if constexpr (contains_waters) {
                        calculator.enqueue_calculate_cross(body1_waters, body2_sym_atomic, 1, cross_merge_id_aw);
                    }
                }
            }
        }
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

template std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate<true>(const data::Molecule&);
template std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate<false>(const data::Molecule&);
template std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate<true,  true>(const data::Molecule&);
template std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate<true,  false>(const data::Molecule&);
template std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate<false, true>(const data::Molecule&);
template std::unique_ptr<hist::ICompositeDistanceHistogram> symmetry::SymmetryManagerMT::calculate<false, false>(const data::Molecule&);