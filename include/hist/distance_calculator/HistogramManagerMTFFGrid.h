#pragma once

#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>

namespace hist {
    template<bool use_weighted_distribution> 
    class HistogramManagerMTFFGrid : public HistogramManagerMTFFAvg<use_weighted_distribution> {
        public:
            using HistogramManagerMTFFAvg<use_weighted_distribution>::HistogramManagerMTFFAvg;

            virtual ~HistogramManagerMTFFGrid() override;

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