#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>

namespace ausaxs::hist {
    class CompositeDistanceHistogramWAXSiS : public CompositeDistanceHistogramFFGrid {
        public:
            using CompositeDistanceHistogramFFGrid::CompositeDistanceHistogramFFGrid;
            CompositeDistanceHistogramWAXSiS(CompositeDistanceHistogramWAXSiS&&) noexcept;
            CompositeDistanceHistogramWAXSiS& operator=(CompositeDistanceHistogramWAXSiS&&) noexcept;
            ~CompositeDistanceHistogramWAXSiS() override;

            /**
             * @brief Override the cache_refresh_intensity_profiles function to change how the exv factor is fitted. 
             */
            void cache_refresh_intensity_profiles(bool sinqd_changed, bool cw_changed, bool cx_changed) const override;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramWAXSiS>, "CompositeDistanceHistogramWAXSiS should support nothrow move semantics.");
}