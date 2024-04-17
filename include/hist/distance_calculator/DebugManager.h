#pragma once

#include <hist/distance_calculator/HistogramManager.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>

namespace hist {
    class DebugDistanceHistogram : public CompositeDistanceHistogramFFAvg {
        public: 
            DebugDistanceHistogram(observer_ptr<const data::Molecule> protein);
            ~DebugDistanceHistogram() override;

            ScatteringProfile debye_transform() const override;

        private:
            observer_ptr<const data::Molecule> protein;
    };

	/**
	 * @brief A multi-threaded simple distance calculator. 
     *        This class is only intended for testing. Use the PartialHistogramManagerMT class for production.
	 */
	template<bool use_weighted_distribution>
	class DebugManager : public HistogramManager<use_weighted_distribution> {
		public:
			using HistogramManager<use_weighted_distribution>::HistogramManager;

			~DebugManager() override;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			std::unique_ptr<DistanceHistogram> calculate() override;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;
	};
}