#pragma once

namespace hist {
    class DistanceHistogram;
    class CompositeDistanceHistogram;
    class CompositeDistanceHistogramFFAvg;
    class CompositeDistanceHistogramFFExplicit;
    class Histogram;
    class HistogramManager;

    /**
     * @brief A ScatteringProfile is just a (q, I(q)) histogram. 
     */    
    using ScatteringProfile = Histogram;
}