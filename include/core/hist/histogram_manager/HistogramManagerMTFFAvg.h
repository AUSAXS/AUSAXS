// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/HistogramManager.h>
#include <hist/detail/CompactCoordinatesFF.h>

namespace ausaxs::hist {
	/**
	 * @brief A histogram manager which uses an average excluded volume approximation. 
	 *
	 * This is equivalent to the CRYSOL implementation, but with a single average excluded volume for all atoms.
	 * To use unique excluded volumes for each atom, see HistogramManagerMTFFExplicit. 
	 */
	template<bool weighted_bins, bool variable_bin_width>
	class HistogramManagerMTFFAvg : public HistogramManager<weighted_bins, variable_bin_width> {
		public:
			using HistogramManager<weighted_bins, variable_bin_width>::HistogramManager;

			virtual ~HistogramManagerMTFFAvg() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

		protected:
			// data stored for inheritance
			std::unique_ptr<hist::detail::CompactCoordinatesFF<variable_bin_width>> data_a_ptr;
		    std::unique_ptr<hist::detail::CompactCoordinatesFF<variable_bin_width>> data_w_ptr;
	};
}