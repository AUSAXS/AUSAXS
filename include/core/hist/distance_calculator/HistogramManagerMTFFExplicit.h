#pragma once

#include <hist/distance_calculator/HistogramManager.h>

namespace ausaxs::hist {
	class CompositeDistanceHistogram;
	namespace detail {class CompactCoordinatesFF;}

	/**
	 * @brief A histogram manager using explicit excluded volume form factors for each atomic type.
	 *		  This is equivalent to the CRYSOL implementation. 
	 */
	template<bool use_weighted_distribution>
	class HistogramManagerMTFFExplicit : public HistogramManager<use_weighted_distribution> {
		public:
			using HistogramManager<use_weighted_distribution>::HistogramManager;

			virtual ~HistogramManagerMTFFExplicit() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

		protected:
			std::unique_ptr<hist::detail::CompactCoordinatesFF> data_a_ptr;
		    std::unique_ptr<hist::detail::CompactCoordinatesFF> data_w_ptr;
	};
}