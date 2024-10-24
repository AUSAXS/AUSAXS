#pragma once

#include <hist/distance_calculator/HistogramManager.h>
#include <hist/detail/CompactCoordinatesFF.h>

namespace ausaxs::hist {
	/**
	 * @brief A histogram manager which uses an average excluded volume approximation. 
	 *
	 * This is equivalent to the CRYSOL implementation, but with a single average excluded volume for all atoms.
	 * To use unique excluded volumes for each atom, see HistogramManagerMTFFExplicit. 
	 */
	template<bool use_weighted_distribution>
	class HistogramManagerMTFFAvg : public HistogramManager<use_weighted_distribution> {
		public:
			using HistogramManager<use_weighted_distribution>::HistogramManager;

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
			std::unique_ptr<hist::detail::CompactCoordinatesFF> data_a_ptr;
		    std::unique_ptr<hist::detail::CompactCoordinatesFF> data_w_ptr;
	};
}