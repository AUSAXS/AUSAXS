// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/HistogramManagerMTBase.h>

namespace ausaxs::hist {
	/**
	 * @brief A histogram manager which uses an average excluded volume approximation. 
	 *
	 * This is equivalent to the CRYSOL implementation, but with a single average excluded volume for all atoms.
	 * To use unique excluded volumes for each atom, see HistogramManagerMTFFExplicit. 
	 */
	template<bool weighted_bins>
	// NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
	class HistogramManagerMTFFAvg : public HistogramManagerMTBase<weighted_bins, true> {
		public:
			using HistogramManagerMTBase<weighted_bins, true>::HistogramManagerMTBase;

			~HistogramManagerMTFFAvg() override;

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
