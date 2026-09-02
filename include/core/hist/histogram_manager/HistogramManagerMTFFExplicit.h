// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/HistFwd.h>
#include <hist/histogram_manager/HistogramManagerMTFFBase.h>

namespace ausaxs::hist {
	/**
	 * @brief A histogram manager using explicit excluded volume form factors for each atomic type.
	 *		  This is equivalent to the CRYSOL implementation. 
	 */
	template<bool weighted_bins, bool variable_bin_width>
	class HistogramManagerMTFFExplicit : public HistogramManagerMTFFBase<weighted_bins, variable_bin_width> {
		public:
			using HistogramManagerMTFFBase<weighted_bins, variable_bin_width>::HistogramManagerMTFFBase;

			virtual ~HistogramManagerMTFFExplicit() override;

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
