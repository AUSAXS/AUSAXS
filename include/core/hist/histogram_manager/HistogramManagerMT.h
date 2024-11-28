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
	template<bool use_weighted_distribution>
	class HistogramManagerMT : public HistogramManager<use_weighted_distribution> {
		public:
			using HistogramManager<use_weighted_distribution>::HistogramManager;
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