#pragma once

#include <hist/distance_calculator/HistogramManagerMT.h>

namespace hist {
    class HistogramManagerMTFFGrid : public HistogramManagerMT {
        public:
            using HistogramManagerMT::HistogramManagerMT;

            ~HistogramManagerMTFFGrid() override;

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