#pragma once

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>

namespace ausaxs::hist {
    /**
     * @brief A histogram manager which uses a grid-based approximation of the excluded volume.
     *        Due to the highly ordered grid structure, weighted bins is required to use this class. 
     */
    class HistogramManagerMTFFGrid : public HistogramManagerMTFFAvg<true> {
        public:
            using HistogramManagerMTFFAvg::HistogramManagerMTFFAvg;

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