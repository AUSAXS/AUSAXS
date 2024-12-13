#include "utility/Exceptions.h"
#include <data/SymmetryManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

#include <cassert>

using namespace ausaxs;
using namespace hist::detail;

enum type {
    SELF_AA, SELF_AW, SELF_WW, 
    SELF_SYM_AA, SELF_SYM_AW, SELF_SYM_WW,
    CROSS_AA, CROSS_AW, CROSS_WW
};

struct _data {
    std::vector<std::vector<CompactCoordinates>> atomic;
    std::vector<std::vector<CompactCoordinates>> waters;
};

struct ScaleResult {
    ScaleResult() = default;
    ScaleResult(int index, int scale) : index(index), scale(scale) {}
    ScaleResult(int index, std::size_t scale) : index(index), scale(static_cast<int>(scale)) {}
    int index;
    int scale;
};

std::vector<_data> generate_transformed_data(const data::Molecule& protein) {
    std::vector<_data> res(protein.size_body());

    // for every body
    for (int i_body1 = 0; i_body1 < static_cast<int>(protein.size_body()); ++i_body1) {
        const auto& body = protein.get_body(i_body1);
        CompactCoordinates data_a(body.get_atoms());
        CompactCoordinates data_w(body.get_waters());

        Vector3<float> cm;
        {
            auto tmp = body.get_cm();
            cm = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]), static_cast<float>(tmp[2])};
        }

        std::vector<std::vector<CompactCoordinates>> atomic(1+body.size_symmetry());
        std::vector<std::vector<CompactCoordinates>> water(1+body.size_symmetry());

        for (int i = 0; i < static_cast<int>(data_a.get_data().size()); ++i) {
            std::cout << "before: " << data_a.get_data()[i].value.pos << std::endl;
        }

        // loop over its symmetries
        for (int i_sym_1 = 0; i_sym_1 < static_cast<int>(body.size_symmetry()); ++i_sym_1) {
            const auto& symmetry = body.get_symmetry(i_sym_1);

            std::vector<CompactCoordinates> sym_atomic(symmetry.repeat, data_a);
            std::vector<CompactCoordinates> sym_water(symmetry.repeat, data_w);

            // for every symmetry, loop over how many times it should be repeated
            // it is then repeatedly applied to the same data
            for (int i_repeat = 0; i_repeat < symmetry.repeat; ++i_repeat) {
                auto t = symmetry.get_transform(cm, i_repeat+1);
                std::transform(
                    sym_atomic[i_repeat].get_data().begin(), 
                    sym_atomic[i_repeat].get_data().end(), 
                    sym_atomic[i_repeat].get_data().begin(), 
                    [t] (const CompactCoordinatesData& v) -> CompactCoordinatesData {return {t(v.value.pos), v.value.w}; }
                );
                std::transform(
                    sym_water[i_repeat].get_data().begin(), 
                    sym_water[i_repeat].get_data().end(), 
                    sym_water[i_repeat].get_data().begin(), 
                    [t] (const CompactCoordinatesData& v) -> CompactCoordinatesData {return {t(v.value.pos), v.value.w}; }
                );
                for (int i = 0; i < static_cast<int>(sym_atomic[i_repeat].get_data().size()); ++i) {
                    std::cout << "after: " << sym_atomic[i_repeat].get_data()[i].value.pos << std::endl;
                }
            }
            atomic[1+i_sym_1] = std::move(sym_atomic);
            water [1+i_sym_1] = std::move(sym_water);
        }
        atomic[0] = {std::move(data_a)};
        water[0] = {std::move(data_w)};
        res[i_body1] = {std::move(atomic), std::move(water)};
    }
    return res;
}

template<bool use_weighted_distribution>
std::unique_ptr<hist::CompositeDistanceHistogram> hist::detail::SymmetryManager::calculate(
    const data::Molecule& protein
) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    hist::distance_calculator::SimpleCalculator<use_weighted_distribution> calculator;

    // start by generating the transformed data
    // note that we are responsible for guaranteeing their lifetime until all enqueue_calculate_* calls are done
    auto data = generate_transformed_data(protein);

    // since the results will be mixed together into a single anonymous vector, we need to keep track of which index corresponds to which type
    // we need this to distinguish between e.g. self and cross histograms
    std::vector<type> self_indices, cross_indices;

    // some of the calculations can be reused multiple times, so we store them in a vector so we can scale them later
    std::vector<ScaleResult> scale_result;
    for (int i_body1 = 0; i_body1 < static_cast<int>(protein.size_body()); ++i_body1) {
        const auto& body = protein.get_body(i_body1);
        const auto& body1_atomic = data[i_body1].atomic[0][0];
        const auto& body1_waters = data[i_body1].waters[0][0];

        // all self calculations must be scaled, so we don't need to store them
        calculator.enqueue_calculate_self(body1_atomic);
        calculator.enqueue_calculate_self(body1_waters); 
        scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_atomic, body1_waters), 1 + body.size_symmetry_total());
        self_indices.emplace_back(SELF_AA);
        self_indices.emplace_back(SELF_WW);
        cross_indices.emplace_back(SELF_AW);

        for (int i_sym1 = 0; i_sym1 < static_cast<int>(body.size_symmetry()); ++i_sym1) {
            const auto& sym1 = body.get_symmetry(i_sym1);
            for (int i_repeat1 = 0; i_repeat1 < sym1.repeat; ++i_repeat1) {
                const auto& body1_sym_atomic = data[i_body1].atomic[1+i_sym1][i_repeat1];
                const auto& body1_sym_waters = data[i_body1].waters[1+i_sym1][i_repeat1];

                // assume we have 3 repeats of symmetry B, so we have the bodies: A B1 B2 B3. Then
                // AB1 == AB2 == AB3, (scale AB1 by 3)
                // AB2 == B1B3        (scale B1B3 by 2)
                // AB3                (scale AB3 by 1)
                int scale = sym1.repeat - i_repeat1;
                scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_atomic, body1_sym_atomic), scale);
                scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_atomic, body1_sym_waters), scale);
                cross_indices.emplace_back(SELF_SYM_AA);
                cross_indices.emplace_back(SELF_SYM_AW);

                scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_waters, body1_sym_waters), scale);
                scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_waters, body1_sym_atomic), scale);
                cross_indices.emplace_back(SELF_SYM_WW);
                cross_indices.emplace_back(SELF_SYM_AW);

                // external histograms with other bodies
                for (int j_body1 = i_body1+1; j_body1 < static_cast<int>(protein.size_body()); ++j_body1) {
                    const auto& body2 = protein.get_body(j_body1);
                    const auto& body2_atomic = data[j_body1].atomic[0][0];
                    const auto& body2_waters = data[j_body1].waters[0][0];

                    calculator.enqueue_calculate_cross(body2_atomic, body1_sym_atomic);
                    calculator.enqueue_calculate_cross(body2_atomic, body1_sym_waters);
                    cross_indices.emplace_back(CROSS_AA);
                    cross_indices.emplace_back(CROSS_AW);

                    calculator.enqueue_calculate_cross(body2_waters, body1_sym_waters);
                    calculator.enqueue_calculate_cross(body2_waters, body1_sym_atomic);
                    cross_indices.emplace_back(CROSS_WW);
                    cross_indices.emplace_back(CROSS_AW);

                    // external histograms with other symmetries in same body
                    for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                        const auto& sym2 = body2.get_symmetry(j_sym1);
                        for (int j_repeat1 = 0; j_repeat1 < sym2.repeat; ++j_repeat1) {
                            const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                            const auto& body2_sym_waters = data[j_body1].waters[1+j_sym1][j_repeat1];

                            calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic);
                            calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_waters);
                            cross_indices.emplace_back(SELF_SYM_AA);
                            cross_indices.emplace_back(SELF_SYM_AW);

                            calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_waters);
                            calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_atomic);
                            cross_indices.emplace_back(SELF_SYM_WW);
                            cross_indices.emplace_back(SELF_SYM_AW);
                        }
                    }
                }

                // internal histogram with other symmetries in same body
                for (int i_sym2 = i_sym1+1; i_sym2 < static_cast<int>(body.size_symmetry()); ++i_sym2) {
                    const auto& sym2 = body.get_symmetry(i_sym2);
                    for (int i_repeat2 = 0; i_repeat2 < sym2.repeat; ++i_repeat2) {
                        const auto& body2_sym_atomic = data[i_body1].atomic[1+i_sym2][i_repeat2];
                        const auto& body2_sym_waters = data[i_body1].waters[1+i_sym2][i_repeat2];

                        calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic);
                        calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_waters);
                        cross_indices.emplace_back(SELF_SYM_AA);
                        cross_indices.emplace_back(SELF_SYM_AW);

                        calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_waters);
                        calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_atomic);
                        cross_indices.emplace_back(SELF_SYM_WW);
                        cross_indices.emplace_back(SELF_SYM_AW);
                    }
                }
            }
        }

        // external histograms with other bodies
        for (int j_body1 = i_body1+1; j_body1 < static_cast<int>(protein.size_body()); ++j_body1) {
            const auto& body2 = protein.get_body(j_body1);
            const auto& body2_atomic = data[j_body1].atomic[0][0];
            const auto& body2_waters = data[j_body1].waters[0][0];

            calculator.enqueue_calculate_cross(body1_atomic, body2_atomic);
            calculator.enqueue_calculate_cross(body1_atomic, body2_waters);
            cross_indices.emplace_back(CROSS_AA);
            cross_indices.emplace_back(CROSS_AW);

            calculator.enqueue_calculate_cross(body1_waters, body2_waters);
            calculator.enqueue_calculate_cross(body1_waters, body2_atomic);
            cross_indices.emplace_back(CROSS_WW);
            cross_indices.emplace_back(CROSS_AW);

            // external histograms with other symmetries in same body
            for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                const auto& sym2 = body2.get_symmetry(j_sym1);
                for (int j_repeat1 = 0; j_repeat1 < sym2.repeat; ++j_repeat1) {
                    const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                    const auto& body2_sym_waters = data[j_body1].waters[1+j_sym1][j_repeat1];

                    calculator.enqueue_calculate_cross(body1_atomic, body2_sym_atomic);
                    calculator.enqueue_calculate_cross(body1_atomic, body2_sym_waters);
                    cross_indices.emplace_back(CROSS_AA);
                    cross_indices.emplace_back(CROSS_AW);

                    calculator.enqueue_calculate_cross(body1_waters, body2_sym_waters);
                    calculator.enqueue_calculate_cross(body1_waters, body2_sym_atomic);
                    cross_indices.emplace_back(CROSS_WW);
                    cross_indices.emplace_back(CROSS_AW);
                }
            }
        }
    }

    auto res = calculator.run();

    assert(res.self.size() == self_indices.size() && "SymmetryManager::calculate: self size mismatch");
    assert(res.cross.size() == cross_indices.size() && "SymmetryManager::calculate: cross size mismatch");
    assert(2*protein.size_body() == res.self.size() && "SymmetryManager::calculate: self size mismatch");

    auto scale_hist = [] (GenericDistribution1D_t& hist, int scale) {
        assert(0 < scale && "SymmetryManager::calculate: scale must be positive");
        if (scale == 1) {return;}
        for (int i = 0; i < static_cast<int>(hist.size()); ++i) {
            hist.set_content(i, hist.get_content(i)*scale);
        }
    };

    // scale all self histograms
    for (int i = 0; i < static_cast<int>(protein.size_body()); ++i) {
        int duplicates = 1 + protein.get_body(i).size_symmetry_total();
        scale_hist(res.self[2*i],   duplicates);
        scale_hist(res.self[2*i+1], duplicates);
    }

    // scale the selected cross histograms
    for (int i = 0; i < static_cast<int>(scale_result.size()); ++i) {
        scale_hist(res.cross[scale_result[i].index], scale_result[i].scale);
    }

    GenericDistribution1D_t p_aa = res.self[0];
    GenericDistribution1D_t p_ww = res.self[1];
    GenericDistribution1D_t p_aw = res.cross[0];

    // add self terms
    for (int i = 1; i < static_cast<int>(protein.size_body()); ++i) {
        p_aa += res.self[2*i];
        p_ww += res.self[2*i+1];
    }

    // add cross terms
    for (int i = 1; i < static_cast<int>(cross_indices.size()); ++i) {
        switch (cross_indices[i]) {
            case SELF_AA:
            case CROSS_AA:
            case SELF_SYM_AA: p_aa += res.cross[i]; break;

            case SELF_WW:
            case CROSS_WW:
            case SELF_SYM_WW: p_ww += res.cross[i]; break;

            case SELF_AW:
            case CROSS_AW:
            case SELF_SYM_AW: p_aw += res.cross[i]; break;
        }
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
        return std::make_unique<CompositeDistanceHistogram>(
            Distribution1D(std::move(p_aa)), 
            Distribution1D(std::move(p_aw)), 
            Distribution1D(std::move(p_ww)), 
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

template std::unique_ptr<hist::CompositeDistanceHistogram> hist::detail::SymmetryManager::calculate<true>(const data::Molecule&);
template std::unique_ptr<hist::CompositeDistanceHistogram> hist::detail::SymmetryManager::calculate<false>(const data::Molecule&);