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

enum h_type {
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

#define DEBUG_MODE false
#if DEBUG_MODE
    #include <iostream>
    #include <iomanip>
#endif
std::vector<_data> generate_transformed_data(const data::Molecule& protein) {
    std::vector<_data> res(protein.size_body());

    // for every body
    for (int i_body1 = 0; i_body1 < static_cast<int>(protein.size_body()); ++i_body1) {
        const auto& body = protein.get_body(i_body1);
        CompactCoordinates data_a(body.get_atoms());
        CompactCoordinates data_w;
        if (body.size_water() != 0) {
            data_w = CompactCoordinates(body.get_waters());
        } 
        hist::detail::SimpleExvModel::apply_simple_excluded_volume(data_a, &protein);

        Vector3<float> cm;
        {
            auto tmp = body.get_cm();
            cm = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]), static_cast<float>(tmp[2])};
        }

        std::vector<std::vector<CompactCoordinates>> atomic(1+body.size_symmetry());
        std::vector<std::vector<CompactCoordinates>> water(1+body.size_symmetry());

        #if DEBUG_MODE
            for (int i = 0; i < static_cast<int>(data_a.get_data().size()); ++i) {
                std::cout << "before: " << data_a.get_data()[i].value.pos << std::endl;
            }
        #endif

        // loop over its symmetries
        for (int i_sym_1 = 0; i_sym_1 < static_cast<int>(body.size_symmetry()); ++i_sym_1) {
            const auto& symmetry = body.symmetry().get(i_sym_1);

            std::vector<CompactCoordinates> sym_atomic(symmetry.repeat, data_a);
            std::vector<CompactCoordinates> sym_water(symmetry.repeat, data_w);

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
                std::transform(
                    sym_water[i_repeat].get_data().begin(), 
                    sym_water[i_repeat].get_data().end(), 
                    sym_water[i_repeat].get_data().begin(), 
                    [t] (const CompactCoordinatesData& v) -> CompactCoordinatesData {return {t(v.value.pos), v.value.w}; }
                );

                #if DEBUG_MODE
                    for (int i = 0; i < static_cast<int>(sym_atomic[i_repeat].get_data().size()); ++i) {
                        std::cout << "after: " << sym_atomic[i_repeat].get_data()[i].value.pos << std::endl;
                    }
                #endif
            }
            atomic[1+i_sym_1] = std::move(sym_atomic);
            water [1+i_sym_1] = std::move(sym_water);
        }
        atomic[0]    = {std::move(data_a)};
        water[0]     = {std::move(data_w)};
        res[i_body1] = {std::move(atomic), std::move(water)};
    }
    return res;
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
    auto data = generate_transformed_data(protein);

    // since the results will be mixed together into a single anonymous vector, we need to keep track of which index corresponds to which type
    // we need this to distinguish between e.g. self and cross histograms
    std::vector<h_type> self_indices, cross_indices;

    // some of the calculations can be reused multiple times, so we store them in a vector so we can scale them later
    std::vector<ScaleResult> scale_result;
    #if DEBUG_MODE
        std::vector<int> cross_scaling_factor;
        std::vector<int> self_scaling_factor;
        std::vector<std::string> self_msg;
        std::vector<int> self_msg_length;
        std::vector<std::string> cross_msg;
        std::vector<int> cross_msg_length;
    #endif
    for (int i_body1 = 0; i_body1 < static_cast<int>(protein.size_body()); ++i_body1) {
        const auto& body = protein.get_body(i_body1);
        const auto& body1_atomic = data[i_body1].atomic[0][0];
        const auto& body1_waters = data[i_body1].waters[0][0];

        // all self calculations must be scaled, so we don't need to store them
        calculator.enqueue_calculate_self(body1_atomic);
        self_indices.emplace_back(SELF_AA);
        if constexpr (contains_waters) {
            calculator.enqueue_calculate_self(body1_waters);
            scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_atomic, body1_waters), 1 + body.size_symmetry_total());
            self_indices.emplace_back(SELF_WW);
            cross_indices.emplace_back(SELF_AW);
        }

        #if DEBUG_MODE
            self_msg.emplace_back("B" + std::to_string(i_body1) + "00 x B" + std::to_string(i_body1) + "00");
            self_scaling_factor.emplace_back(1 + body.size_symmetry_total());
            if constexpr (contains_waters) {
                cross_msg.emplace_back("B" + std::to_string(i_body1) + "00 x B" + std::to_string(i_body1) + "00");
                cross_msg_length.emplace_back(1);
                self_scaling_factor.emplace_back(1 + body.size_symmetry_total());
                cross_scaling_factor.emplace_back(1 + body.size_symmetry_total());
                self_msg_length.emplace_back(2);
            } else {
                self_msg_length.emplace_back(1);
            }
        #endif

        for (int i_sym1 = 0; i_sym1 < static_cast<int>(body.size_symmetry()); ++i_sym1) {
            const auto& sym1 = body.symmetry().get(i_sym1);
            bool closed = sym1.is_closed();
            for (int i_repeat1 = 0; i_repeat1 < (sym1.repeat - closed); ++i_repeat1) {
                const auto& body1_sym_atomic = data[i_body1].atomic[1+i_sym1][i_repeat1];
                const auto& body1_sym_waters = data[i_body1].waters[1+i_sym1][i_repeat1];

                // assume we have 3 repeats of symmetry B, so we have the bodies: A B1 B2 B3. Then
                // AB1 == AB2 == AB3, (scale AB1 by 3)
                // AB2 == B1B3        (scale B1B3 by 2)
                // AB3                (scale AB3 by 1)
                //
                // if the symmetry is closed, AB3 == AB1
                int scale = sym1.repeat - i_repeat1;
                if (i_repeat1 == 0 && closed) {scale += 1;}
                scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_atomic, body1_sym_atomic), scale);
                cross_indices.emplace_back(SELF_SYM_AA);

                if constexpr (contains_waters) {
                    scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_atomic, body1_sym_waters), scale);
                    cross_indices.emplace_back(SELF_SYM_AW);

                    scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_waters, body1_sym_waters), scale);
                    scale_result.emplace_back(calculator.enqueue_calculate_cross(body1_waters, body1_sym_atomic), scale);
                    cross_indices.emplace_back(SELF_SYM_WW);
                    cross_indices.emplace_back(SELF_SYM_AW);
                }

                #if DEBUG_MODE
                    cross_msg.emplace_back("B" + std::to_string(i_body1) + "00 x B" + std::to_string(i_body1) + std::to_string(1+i_sym1) + std::to_string(1+i_repeat1));
                    if constexpr (contains_waters) {
                        cross_msg_length.emplace_back(4);
                        cross_scaling_factor.emplace_back(scale);                
                        cross_scaling_factor.emplace_back(scale);
                        cross_scaling_factor.emplace_back(scale);                
                        cross_scaling_factor.emplace_back(scale);
                    } else {
                        cross_msg_length.emplace_back(1);
                        cross_scaling_factor.emplace_back(scale);
                    }
                #endif

                // external histograms with other bodies
                for (int j_body1 = i_body1+1; j_body1 < static_cast<int>(protein.size_body()); ++j_body1) {
                    const auto& body2 = protein.get_body(j_body1);
                    const auto& body2_atomic = data[j_body1].atomic[0][0];
                    const auto& body2_waters = data[j_body1].waters[0][0];

                    calculator.enqueue_calculate_cross(body2_atomic, body1_sym_atomic);
                    cross_indices.emplace_back(CROSS_AA);

                    if constexpr (contains_waters) {
                        calculator.enqueue_calculate_cross(body2_atomic, body1_sym_waters);
                        cross_indices.emplace_back(CROSS_AW);

                        calculator.enqueue_calculate_cross(body2_waters, body1_sym_waters);
                        calculator.enqueue_calculate_cross(body2_waters, body1_sym_atomic);
                        cross_indices.emplace_back(CROSS_WW);
                        cross_indices.emplace_back(CROSS_AW);
                    }

                    #if DEBUG_MODE
                        cross_msg.emplace_back("B" + std::to_string(j_body1) + "00 x B" + std::to_string(i_body1) + std::to_string(1+i_sym1) + std::to_string(1+i_repeat1));
                        if constexpr (contains_waters) {
                            cross_msg_length.emplace_back(4);
                            cross_scaling_factor.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                        } else {
                            cross_msg_length.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                        }
                    #endif

                    // external histograms with other symmetries in same body
                    for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                        const auto& sym2 = body2.symmetry().get(j_sym1);
                        for (int j_repeat1 = 0; j_repeat1 < sym2.repeat; ++j_repeat1) {
                            const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                            const auto& body2_sym_waters = data[j_body1].waters[1+j_sym1][j_repeat1];

                            calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic);
                            cross_indices.emplace_back(SELF_SYM_AA);

                            if constexpr (contains_waters) {
                                calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_waters);
                                cross_indices.emplace_back(SELF_SYM_AW);

                                calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_waters);
                                calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_atomic);
                                cross_indices.emplace_back(SELF_SYM_WW);
                                cross_indices.emplace_back(SELF_SYM_AW);
                            }

                            #if DEBUG_MODE
                                cross_msg.emplace_back("B" + std::to_string(i_body1) + std::to_string(1+i_sym1) + std::to_string(1+i_repeat1) + " x B" + std::to_string(j_body1) + std::to_string(1+j_sym1) + std::to_string(1+j_repeat1));
                                if constexpr (contains_waters) {
                                    cross_msg_length.emplace_back(4);
                                    cross_scaling_factor.emplace_back(1);
                                    cross_scaling_factor.emplace_back(1);
                                    cross_scaling_factor.emplace_back(1);
                                    cross_scaling_factor.emplace_back(1);
                                } else {
                                    cross_msg_length.emplace_back(1);
                                    cross_scaling_factor.emplace_back(1);
                                }
                            #endif
                        }
                    }
                }

                // internal histogram with other symmetries in same body
                for (int i_sym2 = i_sym1+1; i_sym2 < static_cast<int>(body.size_symmetry()); ++i_sym2) {
                    const auto& sym2 = body.symmetry().get(i_sym2);
                    for (int i_repeat2 = 0; i_repeat2 < sym2.repeat; ++i_repeat2) {
                        const auto& body2_sym_atomic = data[i_body1].atomic[1+i_sym2][i_repeat2];
                        const auto& body2_sym_waters = data[i_body1].waters[1+i_sym2][i_repeat2];

                        calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic);
                        cross_indices.emplace_back(SELF_SYM_AA);

                        if constexpr (contains_waters) {
                            calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_waters);
                            cross_indices.emplace_back(SELF_SYM_AW);

                            calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_waters);
                            calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_atomic);
                            cross_indices.emplace_back(SELF_SYM_WW);
                            cross_indices.emplace_back(SELF_SYM_AW);
                        }

                        #if DEBUG_MODE
                            cross_msg.emplace_back("B" + std::to_string(i_body1) + std::to_string(1+i_sym1) + std::to_string(1+i_repeat1) + " x B" + std::to_string(i_body1) + std::to_string(1+i_sym2) + std::to_string(1+i_repeat2));
                            if constexpr (contains_waters) {
                                cross_msg_length.emplace_back(4);
                                cross_scaling_factor.emplace_back(1);
                                cross_scaling_factor.emplace_back(1);
                                cross_scaling_factor.emplace_back(1);
                                cross_scaling_factor.emplace_back(1);
                            } else {
                                cross_msg_length.emplace_back(1);
                                cross_scaling_factor.emplace_back(1);
                            }
                        #endif
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
            cross_indices.emplace_back(CROSS_AA);

            if constexpr (contains_waters) {
                calculator.enqueue_calculate_cross(body1_atomic, body2_waters);
                cross_indices.emplace_back(CROSS_AW);

                calculator.enqueue_calculate_cross(body1_waters, body2_waters);
                calculator.enqueue_calculate_cross(body1_waters, body2_atomic);
                cross_indices.emplace_back(CROSS_WW);
                cross_indices.emplace_back(CROSS_AW);
            }

            #if DEBUG_MODE
                cross_msg.emplace_back("B" + std::to_string(i_body1) + "00 x B" + std::to_string(j_body1) + "00");
                if constexpr (contains_waters) {
                    cross_msg_length.emplace_back(4);
                    cross_scaling_factor.emplace_back(1);
                    cross_scaling_factor.emplace_back(1);
                    cross_scaling_factor.emplace_back(1);
                    cross_scaling_factor.emplace_back(1);
                } else {
                    cross_msg_length.emplace_back(1);
                    cross_scaling_factor.emplace_back(1);
                }
            #endif

            // external histograms with other symmetries in same body
            for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                const auto& sym2 = body2.symmetry().get(j_sym1);
                for (int j_repeat1 = 0; j_repeat1 < sym2.repeat; ++j_repeat1) {
                    const auto& body2_sym_atomic = data[j_body1].atomic[1+j_sym1][j_repeat1];
                    const auto& body2_sym_waters = data[j_body1].waters[1+j_sym1][j_repeat1];

                    calculator.enqueue_calculate_cross(body1_atomic, body2_sym_atomic);
                    cross_indices.emplace_back(CROSS_AA);

                    if constexpr (contains_waters) {
                        calculator.enqueue_calculate_cross(body1_atomic, body2_sym_waters);
                        cross_indices.emplace_back(CROSS_AW);

                        calculator.enqueue_calculate_cross(body1_waters, body2_sym_waters);
                        calculator.enqueue_calculate_cross(body1_waters, body2_sym_atomic);
                        cross_indices.emplace_back(CROSS_WW);
                        cross_indices.emplace_back(CROSS_AW);
                    }

                    #if DEBUG_MODE
                        cross_msg.emplace_back("B" + std::to_string(i_body1) + "00 x B" + std::to_string(j_body1) + std::to_string(1+j_sym1) + std::to_string(1+j_repeat1));
                        if constexpr (contains_waters) {
                            cross_msg_length.emplace_back(4);
                            cross_scaling_factor.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                        } else {
                            cross_msg_length.emplace_back(1);
                            cross_scaling_factor.emplace_back(1);
                        }
                    #endif
                }
            }
        }
    }

    auto res = calculator.run();

    #if DEBUG_MODE
        assert(res.self.size() == self_scaling_factor.size());
        assert(res.cross.size() == cross_scaling_factor.size());

        std::cout << "\nself results: " << std::endl;
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << constants::axes::d_axis.get_bin_value(i) << " ";
        }
        std::cout << std::endl;
        int msg_idx = 0;
        int msg_shift_idx = self_msg_length[0];
        std::cout << self_msg[0] << std::endl;
        for (int j = 0; j < static_cast<int>(res.self.size()); ++j) {
            if (j == msg_shift_idx) {
                std::cout << self_msg[++msg_idx] << std::endl;
                msg_shift_idx += self_msg_length[msg_idx];
            }
            auto& r = res.self[j];
            for (int i = 0; i < 20; ++i) {
                std::cout << std::setw(4) << self_scaling_factor[j]*r.get_content(i) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "\ncross results: " << std::endl;
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << constants::axes::d_axis.get_bin_value(i) << " ";
        }
        std::cout << std::endl;
        msg_idx = 0;
        msg_shift_idx = cross_msg_length[0];
        std::cout << cross_msg[0] << std::endl;
        for (int j = 0; j < static_cast<int>(res.cross.size()); ++j) {
            if (j == msg_shift_idx) {
                std::cout << cross_msg[++msg_idx] << std::endl;
                msg_shift_idx += cross_msg_length[msg_idx];
            }
            auto& r = res.cross[j];
            for (int i = 0; i < 20; ++i) {
                std::cout << std::setw(4) << cross_scaling_factor[j]*r.get_content(i) << " ";
            }
            std::cout << std::endl;
        }
    #endif

    assert(res.self.size() == self_indices.size() && "SymmetryManager::calculate: self size mismatch");
    assert(res.cross.size() == cross_indices.size() && "SymmetryManager::calculate: cross size mismatch");
    assert((contains_waters ? 2 : 1)*protein.size_body() == res.self.size() && "SymmetryManager::calculate: self size mismatch");

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
        if constexpr (contains_waters) {
            scale_hist(res.self[2*i],   duplicates);
            scale_hist(res.self[2*i+1], duplicates);
        } else {
            scale_hist(res.self[i], duplicates);
        }
    }

    // scale the selected cross histograms
    for (int i = 0; i < static_cast<int>(scale_result.size()); ++i) {
        scale_hist(res.cross[scale_result[i].index], scale_result[i].scale);
    }

    GenericDistribution1D_t p_aa = res.self[0];
    GenericDistribution1D_t p_ww, p_aw;
    if constexpr (contains_waters) {
        p_ww = res.self[1];
        p_aw = res.cross[0];
    } else {
        p_ww = GenericDistribution1D_t(p_aa.size());
        p_aw = GenericDistribution1D_t(p_aa.size());
    }

    // add self terms
    for (int i = 1; i < static_cast<int>(protein.size_body()); ++i) {
        if constexpr (contains_waters) {
            p_aa += res.self[2*i];
            p_ww += res.self[2*i+1];
        } else {
            p_aa += res.self[i];
        }
    }

    // add cross terms
    for (int i = contains_waters ? 1 : 0; i < static_cast<int>(cross_indices.size()); ++i) {
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