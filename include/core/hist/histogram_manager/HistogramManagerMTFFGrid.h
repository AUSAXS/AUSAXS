// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <grid/detail/GridExcludedVolume.h>

namespace ausaxs::hist {
    /**
     * @brief A histogram manager which uses a grid-based approximation of the excluded volume.
     *        Due to the highly ordered grid structure, weighted bins is required to use this class. 
     */
    template<bool variable_bin_width>
    class HistogramManagerMTFFGrid : public HistogramManagerMTFFAvg<true, variable_bin_width> {
        public:
            using HistogramManagerMTFFAvg<true, variable_bin_width>::HistogramManagerMTFFAvg;

            virtual ~HistogramManagerMTFFGrid() override;

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
    };
}