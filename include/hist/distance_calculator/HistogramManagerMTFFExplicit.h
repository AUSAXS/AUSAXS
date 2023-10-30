#pragma once

#include <hist/distance_calculator/HistogramManager.h>

namespace hist {
	class CompositeDistanceHistogram;
	namespace detail {class CompactCoordinatesFF;}

	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
	 */
	class HistogramManagerMTFFExplicit : public HistogramManager {
		public:
			using HistogramManager::HistogramManager;

			HistogramManagerMTFFExplicit(HistogramManager&);

			virtual ~HistogramManagerMTFFExplicit() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<CompositeDistanceHistogram> calculate_all() override;

		protected:
			std::unique_ptr<hist::detail::CompactCoordinatesFF> data_p_ptr;
		    std::unique_ptr<hist::detail::CompactCoordinatesFF> data_h_ptr;
	};
}