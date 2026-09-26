// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <hist/histogram_manager/HistogramManager.h>

#include <memory>
#include <type_traits>
#include <vector>

namespace ausaxs::hist {
	/**
	 * @brief Common machinery for the multithreaded histogram managers which calculate the whole histogram in one go.
	 *        The derived managers only decide what is built from the distance distributions.
	 *
	 * @tparam form_factors Whether the atoms are resolved by form factor. If so, they are split into one set per active form factor
	 *         type and their pairs are only counted, since the form factors are applied later by the intensity calculator.
	 *         Otherwise all atoms are a single set, and each pair is weighted by the product of the weights of its atoms.
	 */
	template<bool weighted_bins, bool form_factors>
	// NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
	class HistogramManagerMTBase : public HistogramManager<weighted_bins> {
		public:
			using HistogramManager<weighted_bins>::HistogramManager;

			~HistogramManagerMTBase() override;

		protected:
			/**
			 * @brief The pairwise distance distributions, trimmed to the bins holding anything, before any excluded volume accounting.
			 *        With form_factors, the atomic distributions are resolved by form factor type.
			 */
			struct Distributions {
				std::conditional_t<form_factors,
					typename GenericDistribution3D<weighted_bins, Shape::Triangular>::type, // unordered (ff_type1, ff_type2), distance
					typename GenericDistribution1D<weighted_bins>::type
				> p_aa;
				std::conditional_t<form_factors,
					typename GenericDistribution2D<weighted_bins>::type,                     // ff_type, distance
					typename GenericDistribution1D<weighted_bins>::type
				> p_aw;
				typename GenericDistribution1D<weighted_bins>::type p_ww;
				typename GenericDistribution1D<weighted_bins>::type p_tot;
			};

			/**
			 * @brief Build the compact coordinates for the atoms and waters, and evaluate all pairwise distances between them, including the self-correlations.
			 */
			Distributions compute_distributions();

			// data stored for inheritance: the atoms, with form_factors split into one set per active type, and the waters whole.
			// with form_factors their pairs are only counted, so any calculation on them must use a unit_weights Calculator
			std::unique_ptr<std::conditional_t<form_factors, std::vector<hist::detail::CompactCoordinates>, hist::detail::CompactCoordinates>> data_a_ptr;
			std::unique_ptr<hist::detail::CompactCoordinates> data_w_ptr;
	};
}
