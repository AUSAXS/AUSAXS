#pragma once

#include <hist/HistogramManager.h>

namespace hist {
	class CompositeDistanceHistogram;

	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
	 */
	class HistogramManagerMTFF : public HistogramManager {
		public:
			using HistogramManager::HistogramManager;

			HistogramManagerMTFF(HistogramManager&);

			~HistogramManagerMTFF() override;

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