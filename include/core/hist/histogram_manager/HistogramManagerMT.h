// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/HistogramManagerMTBase.h>

namespace ausaxs::hist {
	/**
	 * @brief A multi-threaded simple distance calculator.
	 *
	 * This class does not account for the excluded volume in any way.
	 * To implicitly include it, subtract the average excluded volume charge from each atom.
	 */
	template<bool weighted_bins>
	// NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
	class HistogramManagerMT : public HistogramManagerMTBase<weighted_bins, false> {
		public:
			using HistogramManagerMTBase<weighted_bins, false>::HistogramManagerMTBase;
			~HistogramManagerMT() override;

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