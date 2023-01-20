#pragma once

#include <hist/HistogramManager.h>

namespace hist {
	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
	 */
	class HistogramManagerMT : public HistogramManager {
		public:
			using HistogramManager::HistogramManager;

			HistogramManagerMT(HistogramManager&);

			~HistogramManagerMT() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			Histogram calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			ScatteringHistogram calculate_all() override;
	};
}