// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/HistogramManager.h>
#include <hist/detail/CompactCoordinates.h>

namespace ausaxs::hist {
	/**
	 * @brief A multi-threaded simple distance calculator. 
	 *
	 * This class does not account for the excluded volume in any way. 
	 * To implicitly include it, subtract the average excluded volume charge from each atom. 
	 */
	template<bool weighted_bins, bool variable_bin_width>
	class HistogramManagerMT : public HistogramManager<weighted_bins, variable_bin_width> {
		public:
			using HistogramManager<weighted_bins, variable_bin_width>::HistogramManager;
			virtual ~HistogramManagerMT() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;
	};
}