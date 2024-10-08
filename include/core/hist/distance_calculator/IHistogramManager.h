#pragma once

#include <hist/detail/BodyTracker.h>
#include <hist/HistFwd.h>

#include <memory>

namespace hist {
	/**
	 * @brief A generic histogram manager interface. 
     *        This class is responsible for calculating all distances for a given molecule.
	 */
	class IHistogramManager : public hist::BodyTracker {
		public:
            using BodyTracker::BodyTracker;
			virtual ~IHistogramManager() = default;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual std::unique_ptr<DistanceHistogram> calculate() = 0;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual std::unique_ptr<ICompositeDistanceHistogram> calculate_all() = 0;
    };
}