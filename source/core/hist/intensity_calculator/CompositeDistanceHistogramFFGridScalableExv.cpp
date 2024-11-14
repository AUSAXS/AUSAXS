#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>

using namespace hist;

CompositeDistanceHistogramFFGridScalableExv::CompositeDistanceHistogramFFGridScalableExv(
    CompositeDistanceHistogramFFGrid&& res, 
    std::function<std::unique_ptr<ICompositeDistanceHistogram>(double)> eval_scaled_exv
) : CompositeDistanceHistogramFFGrid(std::move(res)), eval_scaled_exv(eval_scaled_exv) {}