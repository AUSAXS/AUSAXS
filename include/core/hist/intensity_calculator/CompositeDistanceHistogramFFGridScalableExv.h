#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <table/VectorDebyeTable.h>
#include <utility/TypeTraits.h>

namespace hist {
    class CompositeDistanceHistogramFFGridScalableExv : public CompositeDistanceHistogramFFGrid {
        public:
            CompositeDistanceHistogramFFGridScalableExv(CompositeDistanceHistogramFFGrid&& cdh, std::function<std::unique_ptr<CompositeDistanceHistogramFFGrid>(double)> eval_scaled_exv);
            void apply_excluded_volume_scaling_factor(double k) override;
            virtual double exv_factor(double) const override {return 1;}

        private:
            std::function<std::unique_ptr<CompositeDistanceHistogramFFGrid>(double)> eval_scaled_exv;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFGridScalableExv>, "CompositeDistanceHistogramFFGridScalableExv should support nothrow move semantics.");
}