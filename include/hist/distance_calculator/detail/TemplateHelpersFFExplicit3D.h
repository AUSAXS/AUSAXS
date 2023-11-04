#pragma once

#include <hist/distance_calculator/detail/TemplateHelpersFFAvg3D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactorType.h>

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
inline void evaluate8(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_aa, typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_ax, typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_xx, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff3::add8::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    for (unsigned int k = 0; k < 8; ++k) {
        p_aa.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], factor*res.weights[k]);
        p_ax.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], factor*data_j[j+k].value.w);
        p_xx.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], factor*1);
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
inline void evaluate4(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_aa, typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_ax, typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_xx, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff3::add4::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    for (unsigned int k = 0; k < 4; ++k) {
        p_aa.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], factor*res.weights[k]);
        p_ax.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], factor*data_j[j+k].value.w);
        p_xx.add(data_i.get_ff_type(i), data_j.get_ff_type(j+k), res.distances[k], factor*1);
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
inline void evaluate1(typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_aa, typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_ax, typename hist::GenericDistribution3D<use_weighted_distribution>::type& p_xx, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinatesFF& data_j, int i, int j) {
    auto res = detail::ff3::add1::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
    p_aa.add(data_i.get_ff_type(i), data_j.get_ff_type(j), res.distance, factor*res.weight);
    p_ax.add(data_i.get_ff_type(i), data_j.get_ff_type(j), res.distance, factor*data_j[j].value.w);
    p_xx.add(data_i.get_ff_type(i), data_j.get_ff_type(j), res.distance, factor*1);
}