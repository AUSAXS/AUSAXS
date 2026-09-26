// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/HistFwd.h>
#include <hist/histogram_manager/HistogramManagerMTBase.h>
#include <settings/ExvSettings.h>

namespace ausaxs::hist {
	/**
	 * @brief A histogram manager using explicit excluded volume form factors for each atomic type.
	 *		  This is equivalent to the CRYSOL implementation. 
	 */
	template<bool weighted_bins>
	// NOLINTNEXTLINE - the destructor is virtual through the dependent base, which the check cannot see on the template pattern
	class HistogramManagerMTFFExplicit : public HistogramManagerMTBase<weighted_bins, true> {
		public:
			/**
			 * @param exv_method The explicit excluded volume model to build the result for; see detail::make_explicit_histogram.
			 */
			explicit HistogramManagerMTFFExplicit(observer_ptr<const data::Molecule> protein, settings::exv::ExvMethod exv_method = settings::exv::exv_method)
				: HistogramManagerMTBase<weighted_bins, true>(protein), exv_method(exv_method) {}

			~HistogramManagerMTFFExplicit() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

		private:
			settings::exv::ExvMethod exv_method;
	};
}
