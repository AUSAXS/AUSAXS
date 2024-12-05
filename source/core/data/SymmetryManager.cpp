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
    std::vector<CompactCoordinates> atomic;
    std::vector<CompactCoordinates> waters;
};

std::vector<_data> generate_transformed_data(const data::Molecule& protein) {
    std::vector<_data> res(protein.size_body());
    for (unsigned int i = 0; i < protein.size_body(); ++i) {
        const auto& body = protein.get_body(i);
        CompactCoordinates data_a(body.get_atoms());
        CompactCoordinates data_w(body.get_waters());
        auto cm = body.get_cm();

        std::vector<CompactCoordinates> atomic(1+body.size_symmetry());
        std::vector<CompactCoordinates> water(1+body.size_symmetry());
        for (unsigned int j = 0; j < body.size_symmetry(); ++j) {
            const auto& symmetry = body.get_symmetries()[j];
            Matrix<float> rotate = matrix::rotation_matrix<float>(symmetry.rotate.x(), symmetry.rotate.y(), symmetry.rotate.z());
            auto transform = [&symmetry, &rotate, &cm] (const CompactCoordinatesData& a) -> CompactCoordinatesData {
                return {rotate*(a.value.pos-cm) + cm + symmetry.translate, a.value.w};
            };

            CompactCoordinates data_a_transformed(data_a);
            CompactCoordinates data_w_transformed(data_w);
            std::transform(
                data_a_transformed.get_data().begin(), 
                data_a_transformed.get_data().end(), 
                data_a_transformed.get_data().begin(), 
                transform
            );

            std::transform(
                data_w_transformed.get_data().begin(), 
                data_w_transformed.get_data().end(), 
                data_w_transformed.get_data().begin(), 
                transform
            );

            atomic[j+1] = std::move(data_a_transformed);
            water[j+1] = std::move(data_w_transformed);
        }
        atomic[0] = std::move(data_a);
        water[0] = std::move(data_w);
        res[i] = {std::move(atomic), std::move(water)};
    }
    return res;
}

template<bool use_weighted_distribution>
std::unique_ptr<hist::CompositeDistanceHistogram> hist::detail::SymmetryManager::calculate(
    const data::Molecule& protein
) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    hist::distance_calculator::SimpleCalculator<use_weighted_distribution> calculator;

    auto data = generate_transformed_data(protein);
    std::vector<type> self_indices, cross_indices;
    std::vector<int> needs_scaling;
    for (unsigned int i = 0; i < protein.size_body(); ++i) {
        const auto& body = protein.get_body(i);
        const auto& body1_atomic = data[i].atomic[0];
        const auto& body1_waters = data[i].waters[0];

        calculator.enqueue_calculate_self(body1_atomic);
        calculator.enqueue_calculate_self(body1_waters);
        needs_scaling.push_back(calculator.enqueue_calculate_cross(body1_atomic, body1_waters));
        self_indices.push_back(SELF_AA);
        self_indices.push_back(SELF_WW);
        cross_indices.push_back(SELF_AW);

        for (unsigned int j = 0; j < body.size_symmetry(); ++j) {
            const auto& body1_sym_atomic = data[i].atomic[j+1];
            const auto& body1_sym_waters = data[i].waters[j+1];

            calculator.enqueue_calculate_cross(body1_atomic, body1_sym_atomic);
            calculator.enqueue_calculate_cross(body1_atomic, body1_sym_waters);
            cross_indices.push_back(SELF_SYM_AA);
            cross_indices.push_back(SELF_SYM_AW);

            calculator.enqueue_calculate_cross(body1_waters, body1_sym_waters);
            calculator.enqueue_calculate_cross(body1_waters, body1_sym_atomic);
            cross_indices.push_back(SELF_SYM_WW);
            cross_indices.push_back(SELF_SYM_AW);

            // external histograms with other bodies
            for (unsigned int k = i+1; k < protein.size_body(); ++k) {
                const auto& body2 = protein.get_body(k);
                const auto& body2_atomic = data[k].atomic[0];
                const auto& body2_waters = data[k].waters[0];

                calculator.enqueue_calculate_cross(body2_atomic, body1_sym_atomic);
                calculator.enqueue_calculate_cross(body2_atomic, body1_sym_waters);
                cross_indices.push_back(CROSS_AA);
                cross_indices.push_back(CROSS_AW);

                calculator.enqueue_calculate_cross(body2_waters, body1_sym_waters);
                calculator.enqueue_calculate_cross(body2_waters, body1_sym_atomic);
                cross_indices.push_back(CROSS_WW);
                cross_indices.push_back(CROSS_AW);

                // external histograms with other symmetries in same body
                for (unsigned int l = 0; l < body2.size_symmetry(); ++l) {
                    const auto& body2_sym_atomic = data[k].atomic[l+1];
                    const auto& body2_sym_waters = data[k].waters[l+1];

                    calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic);
                    calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_waters);
                    cross_indices.push_back(SELF_SYM_AA);
                    cross_indices.push_back(SELF_SYM_AW);

                    calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_waters);
                    calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_atomic);
                    cross_indices.push_back(SELF_SYM_WW);
                    cross_indices.push_back(SELF_SYM_AW);
                }
            }

            // internal histogram with other symmetries in same body
            for (unsigned int k = j+1; k < body.size_symmetry(); ++k) {
                const auto& body2_sym_atomic = data[i].atomic[k+1];
                const auto& body2_sym_waters = data[i].waters[k+1];

                calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic);
                calculator.enqueue_calculate_cross(body1_sym_atomic, body2_sym_waters);
                cross_indices.push_back(SELF_SYM_AA);
                cross_indices.push_back(SELF_SYM_AW);

                calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_waters);
                calculator.enqueue_calculate_cross(body1_sym_waters, body2_sym_atomic);
                cross_indices.push_back(SELF_SYM_WW);
                cross_indices.push_back(SELF_SYM_AW);
            }
        }

        // external histograms with other bodies
        for (unsigned int j = i+1; j < protein.size_body(); ++j) {
            const auto& body2 = protein.get_body(j);
            const auto& body2_atomic = data[j].atomic[0];
            const auto& body2_waters = data[j].waters[0];

            calculator.enqueue_calculate_cross(body1_atomic, body2_atomic);
            calculator.enqueue_calculate_cross(body1_atomic, body2_waters);
            cross_indices.push_back(CROSS_AA);
            cross_indices.push_back(CROSS_AW);

            calculator.enqueue_calculate_cross(body1_waters, body2_waters);
            calculator.enqueue_calculate_cross(body1_waters, body2_atomic);
            cross_indices.push_back(CROSS_WW);
            cross_indices.push_back(CROSS_AW);

            // external histograms with other symmetries in same body
            for (unsigned int k = 0; k < body2.size_symmetry(); ++k) {
                const auto& body2_sym_atomic = data[j].atomic[k+1];
                const auto& body2_sym_waters = data[j].waters[k+1];

                calculator.enqueue_calculate_cross(body1_atomic, body2_sym_atomic);
                calculator.enqueue_calculate_cross(body1_atomic, body2_sym_waters);
                cross_indices.push_back(CROSS_AA);
                cross_indices.push_back(CROSS_AW);

                calculator.enqueue_calculate_cross(body1_waters, body2_sym_waters);
                calculator.enqueue_calculate_cross(body1_waters, body2_sym_atomic);
                cross_indices.push_back(CROSS_WW);
                cross_indices.push_back(CROSS_AW);
            }
        }
    }

    auto res = calculator.run();

    assert(res.self.size() == self_indices.size() && "SymmetryManager::calculate: self size mismatch");
    assert(res.cross.size() == cross_indices.size() && "SymmetryManager::calculate: cross size mismatch");

    GenericDistribution1D_t p_aa = res.self[0];
    GenericDistribution1D_t p_ww = res.self[1];
    GenericDistribution1D_t p_aw = res.cross[0];
    {
        int duplicates = 1 + protein.get_body(0).size_symmetry();
        for (int i = 0; i < static_cast<int>(p_aa.size()); ++i) {
            p_aa.set_content(i, p_aa.get_content(i)*duplicates);
            p_ww.set_content(i, p_ww.get_content(i)*duplicates);
            p_aw.set_content(i, p_aw.get_content(i)*duplicates);
        }
    }

    for (int i = 1; i < static_cast<int>(protein.size_body()); ++i) {
        int duplicates = 1 + protein.get_body(i).size_symmetry();
        for (int j = 0; j < static_cast<int>(p_aa.size()); ++j) {
            res.self[2*i  ].set_content(j, res.self[2*i  ].get_content(j)*duplicates);
            res.self[2*i+1].set_content(j, res.self[2*i+1].get_content(j)*duplicates);
            res.cross[needs_scaling[i]].set_content(j, res.cross[needs_scaling[i]].get_content(j)*duplicates);
        }
        p_aa += res.self[2*i];
        p_ww += res.self[2*i+1];
    }

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