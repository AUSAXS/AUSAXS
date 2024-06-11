#pragma once

#include <hist/distance_calculator/HistogramManager.h>
#include <hist/detail/CompactCoordinates.h>

namespace hist {
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

		protected:
			std::unique_ptr<hist::detail::CompactCoordinates> data_a_ptr;
		    std::unique_ptr<hist::detail::CompactCoordinates> data_w_ptr;
	};
}