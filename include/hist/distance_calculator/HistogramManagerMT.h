#pragma once

#include <hist/distance_calculator/HistogramManager.h>
#include <hist/detail/CompactCoordinates.h>

namespace hist {
	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
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

		protected:
			std::unique_ptr<hist::detail::CompactCoordinates> data_p_ptr;
		    std::unique_ptr<hist::detail::CompactCoordinates> data_h_ptr;
	};
}