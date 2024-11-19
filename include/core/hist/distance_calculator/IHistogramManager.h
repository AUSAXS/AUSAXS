#pragma once

#include <hist/detail/BodyTracker.h>
#include <hist/HistFwd.h>

#include <memory>

namespace ausaxs::hist {
	/**
	 * @brief A generic histogram manager interface. 
     *        This class is responsible for calculating all distances for a given molecule.
	 */
	class IHistogramManager {
		public:
			virtual ~IHistogramManager() = default;

			/**
			 * @brief Whether to account for the excluded volume by subtracting the average excluded volume charge from each atom.
			 */
			static void use_simple_excluded_volume(bool enabled);

			/**
			 * @brief Signal that the hydration layer was modified. 
			 */
			virtual void signal_modified_hydration_layer() {return;};

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