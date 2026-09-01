// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/HistogramManager.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>

namespace ausaxs::hist {
	/**
	 * @brief Common machinery for the form factor-aware multithreaded histogram managers.
	 */
	template<bool weighted_bins, bool variable_bin_width>
	class HistogramManagerMTFFBase : public HistogramManager<weighted_bins, variable_bin_width> {
		public:
			using HistogramManager<weighted_bins, variable_bin_width>::HistogramManager;

			virtual ~HistogramManagerMTFFBase() override;

		protected:
			/**
			 * @brief The raw pairwise distance distributions, before any excluded volume accounting.
			 */
			struct RawDistributions {
				typename GenericDistribution3D<weighted_bins>::type p_aa; // ff_type1, ff_type2, distance
				typename GenericDistribution2D<weighted_bins>::type p_aw; // ff_type, distance
				typename GenericDistribution1D<weighted_bins>::type p_ww; // distance
				typename GenericDistribution1D<weighted_bins>::type p_tot;
				unsigned int max_bin;
			};

			/**
			 * @brief Build the compact coordinates for the atoms and waters, and evaluate all pairwise distances between them, including the self-correlations.
			 */
			RawDistributions compute_raw_distributions();

			// data stored for inheritance
			std::unique_ptr<hist::detail::CompactCoordinatesFF<variable_bin_width>> data_a_ptr;
		    std::unique_ptr<hist::detail::CompactCoordinatesFF<variable_bin_width>> data_w_ptr;
	};
}
