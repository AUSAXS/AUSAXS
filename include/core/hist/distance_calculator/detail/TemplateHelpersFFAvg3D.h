// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactorType.h>

namespace ausaxs {
    /**
     * @brief Calculate the distances between eight atoms and template add them to the histogram. The contribution of excluded volume dummy atoms is template added to the excluded volume bin.
     * 
     * @tparam use_weighted_distribution Whether to keep track of the distances template added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to template add the distances to.
     * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
     * @param data_j The second atom. The form factor index of this will be used for the second axis of the histogram.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool use_weighted_distribution, int factor>
    inline void evaluate8(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        auto res = detail::add8::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {
                p.template add<factor>(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], res.weights[k]);
                p.template add<factor>(data_i.get_ff_type(i), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distances[k], data_j[j+k].value.w);
                p.template add<factor>(static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distances[k], 1);
        }
    }

    /**
     * @brief Calculate the distances between four atoms and template add them to the histogram. The contribution of excluded volume dummy atoms is template added to the excluded volume bin.
     * 
     * @tparam use_weighted_distribution Whether to keep track of the distances template added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to template add the distances to.
     * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
     * @param data_j The second atom. The form factor index of this will be used for the second axis of the histogram.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool use_weighted_distribution, int factor>
    inline void evaluate4(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        auto res = detail::add4::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {
            p.template add<factor>(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], res.weights[k]);
            p.template add<factor>(data_i.get_ff_type(i), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distances[k], data_j[j+k].value.w);
            p.template add<factor>(static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distances[k], 1);
        }
    }

    /**
     * @brief Calculate the distances between two atoms and template add them to the histogram. The contribution of excluded volume dummy atoms is template added to the excluded volume bin.
     * 
     * @tparam use_weighted_distribution Whether to keep track of the distances template added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to template add the distances to.
     * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
     * @param data_j The second atom. The form factor index of this will be used for the second axis of the histogram.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool use_weighted_distribution, int factor>
    inline void evaluate1(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        auto res = detail::add1::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        p.template add<factor>(data_i.get_ff_type(i), data_j.get_ff_type(j), res.distance, res.weight);
        p.template add<factor>(data_i.get_ff_type(i), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance, data_j[j].value.w);
        p.template add<factor>(static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance, 1);
    }
}