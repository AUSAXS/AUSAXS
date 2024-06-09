#pragma once

#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>

namespace hist {
    class HistogramManagerMTFFGridSurface : public HistogramManagerMTFFAvg<true> {
        public:
            using HistogramManagerMTFFAvg::HistogramManagerMTFFAvg;

            virtual ~HistogramManagerMTFFGridSurface() override;

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