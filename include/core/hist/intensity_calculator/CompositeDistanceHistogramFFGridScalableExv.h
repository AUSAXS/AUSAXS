// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    class CompositeDistanceHistogramFFGridScalableExv : public CompositeDistanceHistogramFFGrid {
        public:
            CompositeDistanceHistogramFFGridScalableExv(CompositeDistanceHistogramFFGrid&& res, std::function<std::unique_ptr<CompositeDistanceHistogramFFGrid>(double)> eval_scaled_exv);
            void apply_excluded_volume_scaling_factor(double k) override;
            double exv_factor(double /*q*/) const override {return 1;}

        private:
            std::function<std::unique_ptr<CompositeDistanceHistogramFFGrid>(double)> eval_scaled_exv;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFGridScalableExv>, "CompositeDistanceHistogramFFGridScalableExv should support nothrow move semantics.");
}