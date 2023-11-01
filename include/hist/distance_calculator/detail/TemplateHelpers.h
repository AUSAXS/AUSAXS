#pragma once

#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>

namespace add8_impl {
    template<int factor>
    inline void add8_impl(hist::Distribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d, const hist::detail::CompactCoordinatesData& e, const hist::detail::CompactCoordinatesData& f, const hist::detail::CompactCoordinatesData& g, const hist::detail::CompactCoordinatesData& h) {
        auto res = data.evaluate_rounded(a, b, c, d, e, f, g, h);
        for (unsigned int k = 0; k < 8; ++k) {p.add(res.distance[k], factor*res.weight[k]);}
    }

    template<int factor>
    inline void add8_impl(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d, const hist::detail::CompactCoordinatesData& e, const hist::detail::CompactCoordinatesData& f, const hist::detail::CompactCoordinatesData& g, const hist::detail::CompactCoordinatesData& h) {
        auto res = data.evaluate(a, b, c, d, e, f, g, h);
        for (unsigned int k = 0; k < 8; ++k) {p.add(res.distance[k], factor*res.weight[k]);}
    }
}

namespace add4_impl {
    template<int factor>
    inline void add4_impl(hist::Distribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d) {
        auto res = data.evaluate_rounded(a, b, c, d);
        for (unsigned int k = 0; k < 4; ++k) {p.add(res.distance[k], factor*res.weight[k]);}
    }

    template<int factor>
    inline void add4_impl(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d) {
        auto res = data.evaluate(a, b, c, d);
        for (unsigned int k = 0; k < 4; ++k) {p.add(res.distance[k], factor*res.weight[k]);}
    }
}

namespace add1_impl {
    template<int factor>
    inline void add1_impl(hist::Distribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a) {
        auto res = data.evaluate_rounded(a);
        p.add(res.distance, factor*res.weight);
    }

    template<int factor>
    inline void add1_impl(hist::WeightedDistribution1D& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a) {
        auto res = data.evaluate(a);
        p.add(res.distance, factor*res.weight);
    }
}

/**
 * @brief Calculate the distances between eight atoms and add them to the histogram.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data The first atom.
 */
template<bool use_weighted_distribution, int factor>
inline void add8(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d, const hist::detail::CompactCoordinatesData& e, const hist::detail::CompactCoordinatesData& f, const hist::detail::CompactCoordinatesData& g, const hist::detail::CompactCoordinatesData& h) {
    add8_impl::add8_impl<factor>(p, data, a, b, c, d, e, f, g, h);
}

/**
 * @brief Calculate the distances between four atoms and add them to the histogram.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data The first atom.
 */
template<bool use_weighted_distribution, int factor>
inline void add4(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a, const hist::detail::CompactCoordinatesData& b, const hist::detail::CompactCoordinatesData& c, const hist::detail::CompactCoordinatesData& d) {
    add4_impl::add4_impl<factor>(p, data, a, b, c, d);
}

/**
 * @brief Calculate the distances between two atoms and add them to the histogram.
 * 
 * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
 * @tparam factor A multiplicative factor for the atomic weights. 
 * @param p The histogram to add the distances to.
 * @param data The first atom.
 */
template<bool use_weighted_distribution, int factor>
inline void add1(typename hist::GenericDistribution1D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesData& data, const hist::detail::CompactCoordinatesData& a) {
    add1_impl::add1_impl<factor>(p, data, a);
}