#pragma once

#include <hist/distance_calculator/HistogramManager.h>
#include <hist/CompositeDistanceHistogramFF.h>

namespace hist {
    class DebugDistanceHistogram : public CompositeDistanceHistogramFF {
        public: 
            DebugDistanceHistogram(const data::Molecule* const protein);
            ~DebugDistanceHistogram() override;

            ScatteringProfile debye_transform() const override;

        private:
            const data::Molecule* const protein;
    };

	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
	 */
	class DebugManager : public HistogramManager {
		public:
			using HistogramManager::HistogramManager;

			~DebugManager() override;

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