#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>

using namespace hist;

CompositeDistanceHistogramFF::~CompositeDistanceHistogramFF() = default;

Histogram CompositeDistanceHistogramFF::debye_transform() const {
    return Histogram();
}

Histogram CompositeDistanceHistogramFF::debye_transform(const std::vector<double>& q) const {
    return Histogram();
}