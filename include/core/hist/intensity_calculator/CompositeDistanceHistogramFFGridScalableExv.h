#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <table/VectorDebyeTable.h>
#include <utility/TypeTraits.h>

namespace hist {
    class CompositeDistanceHistogramFFGridScalableExv : public CompositeDistanceHistogramFFGrid {
        public:
            CompositeDistanceHistogramFFGridScalableExv(CompositeDistanceHistogramFFGrid&& cdh, std::function<std::unique_ptr<ICompositeDistanceHistogram>(double)> eval_scaled_exv);
            std::function<std::unique_ptr<ICompositeDistanceHistogram>(double)> eval_scaled_exv;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFGridScalableExv>, "CompositeDistanceHistogramFFGridScalableExv should support nothrow move semantics.");
}