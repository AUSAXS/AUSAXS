#pragma once

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <grid/detail/GridExcludedVolume.h>

namespace ausaxs::hist {
    /**
     * @brief A histogram manager which uses a grid-based approximation of the excluded volume.
     *        Due to the highly ordered grid structure, weighted bins is required to use this class. 
     */
    class HistogramManagerMTFFGridScalableExv : public HistogramManagerMTFFAvg<true> {
        public:
            using HistogramManagerMTFFAvg::HistogramManagerMTFFAvg;

            virtual ~HistogramManagerMTFFGridScalableExv() override;

            /**
             * @brief Calculate only the total scattering histogram. 
             */
            std::unique_ptr<DistanceHistogram> calculate() override;

            /**
             * @brief Calculate all contributions to the scattering histogram. 
             */
            std::unique_ptr<ICompositeDistanceHistogram> calculate_all() override;

        private:
            virtual grid::exv::GridExcludedVolume get_exv() const;

            double exv_factor = 1;
    };
}