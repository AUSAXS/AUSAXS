#pragma once

#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinatesFF.h>

namespace detail::add8 {
    template<int factor>
    inline auto evaluate(hist::Distribution1D& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
    }

    template<int factor>
    inline auto evaluate(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
    }
}

namespace detail::add4 {
    template<int factor>
    inline auto evaluate(hist::Distribution1D& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
    }

    template<int factor>
    inline auto evaluate(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
    }
}

namespace detail::add1 {
    template<int factor>
    inline auto evaluate(hist::Distribution1D& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
        return data_i[i].evaluate_rounded(data_j[j]);
    }

    template<int factor>
    inline auto evaluate(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
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
inline void add8(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::add8::evaluate<factor>(p, data_i, data_j, i, j);
    for (unsigned int k = 0; k < 8; ++k) {p.add(res.distance[k], factor*res.weight[k]);}
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
inline void add4(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::add4::evaluate<factor>(p, data_i, data_j, i, j);
    for (unsigned int k = 0; k < 4; ++k) {p.add(res.distance[k], factor*res.weight[k]);}
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
inline void add1(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::add1::evaluate<factor>(p, data_i, data_j, i, j);
    p.add(res.distance, factor*res.weight);
}