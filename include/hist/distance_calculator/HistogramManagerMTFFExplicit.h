#pragma once

#include <hist/distance_calculator/HistogramManager.h>

namespace hist {
	class CompositeDistanceHistogram;

	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
	 */
	class HistogramManagerMTFFExplicit : public HistogramManager {
		public:
			using HistogramManager::HistogramManager;

			HistogramManagerMTFFExplicit(HistogramManager&);

			~HistogramManagerMTFFExplicit() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<CompositeDistanceHistogram> calculate_all() override;
	};
}