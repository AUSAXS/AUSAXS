#pragma once

#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinatesFF.h>

namespace detail::ff1::add8 {
    template<int use_weighted_distribution>
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

namespace detail::ff1::add4 {
    template<int use_weighted_distribution>
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

namespace detail::ff1::add1 {
    template<int use_weighted_distribution>
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
 * @brief Calculate the distances between eight atoms and add them to the histogram. No excluded volume bin is added.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data_i The first atom.
 * @param data_j The second atom.
 * @param i The index of the first atom.
 * @param j The index of the second atom.
 */
template<bool use_weighted_distribution, int factor>
inline void evaluate8(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff1::add8::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    for (unsigned int k = 0; k < 8; ++k) {p.add(res.distances[k], factor*res.weights[k]);}
}

/**
 * @brief Calculate the distances between four atoms and add them to the histogram. No excluded volume bin is added.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data_i The first atom.
 * @param data_j The second atom.
 * @param i The index of the first atom.
 * @param j The index of the second atom.
 */
template<bool use_weighted_distribution, int factor>
inline void evaluate4(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff1::add4::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    for (unsigned int k = 0; k < 4; ++k) {p.add(res.distances[k], factor*res.weights[k]);}
}

/**
 * @brief Calculate the distances between two atoms and add them to the histogram. No excluded volume bin is added.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data_i The first atom.
 * @param data_j The second atom.
 * @param i The index of the first atom.
 * @param j The index of the second atom.
 */
template<bool use_weighted_distribution, int factor>
inline void evaluate1(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff1::add1::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    p.add(res.distance, factor*res.weight);
}