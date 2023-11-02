#pragma once

#include <hist/distribution/GenericDistribution3D.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactorType.h>

namespace detail::ff3::add8 {
    template<bool use_weighted_distribution>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j);

    template<>
    inline auto evaluate<false>(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
    }

    template<>
    inline auto evaluate<true>(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
    }
}

namespace detail::ff3::add4 {
    template<bool use_weighted_distribution>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j);

    template<>
    inline auto evaluate<false>(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
    }

    template<>
    inline auto evaluate<true>(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
    }
}

namespace detail::ff3::add1 {
    template<bool use_weighted_distribution>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j);

    template<>
    inline auto evaluate<false>(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate_rounded(data_j[j]);
    }

    template<>
    inline auto evaluate<true>(const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate(data_j[j]);
    }
}

/**
 * @brief Calculate the distances between eight atoms and add them to the histogram. The contribution of excluded volume dummy atoms is added to the excluded volume bin.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
 * @param data_j The second atom. The form factor index of this will be used for the second axis of the histogram.
 * @param i The index of the first atom.
 * @param j The index of the second atom.
 */
template<bool use_weighted_distribution, int factor>
inline void evaluate8(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff3::add8::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    for (unsigned int k = 0; k < 8; ++k) {
        p.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distance[k], factor*res.weight[k]);
        p.add(data_i.get_ff_type(i), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance[k], factor*data_j[j+k].value.w);
        p.add(static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance[k], factor*1);
    }
}

/**
 * @brief Calculate the distances between four atoms and add them to the histogram. The contribution of excluded volume dummy atoms is added to the excluded volume bin.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
 * @param data_j The second atom. The form factor index of this will be used for the second axis of the histogram.
 * @param i The index of the first atom.
 * @param j The index of the second atom.
 */
template<bool use_weighted_distribution, int factor>
inline void evaluate4(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff3::add4::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    for (unsigned int k = 0; k < 4; ++k) {
        p.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distance[k], factor*res.weight[k]);
        p.add(data_i.get_ff_type(i), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance[k], factor*data_j[j+k].value.w);
        p.add(static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance[k], factor*1);
    }
}

/**
 * @brief Calculate the distances between two atoms and add them to the histogram. The contribution of excluded volume dummy atoms is added to the excluded volume bin.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
 * @param data_j The second atom. The form factor index of this will be used for the second axis of the histogram.
 * @param i The index of the first atom.
 * @param j The index of the second atom.
 */
template<bool use_weighted_distribution, int factor>
inline void evaluate1(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff3::add1::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    p.add(data_i.get_ff_type(i), data_j.get_ff_type(j), res.distance, factor*res.weight);
    p.add(data_i.get_ff_type(i), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance, factor*data_j[j].value.w);
    p.add(static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), static_cast<unsigned int>(form_factor::form_factor_t::EXCLUDED_VOLUME), res.distance, factor*1);
}