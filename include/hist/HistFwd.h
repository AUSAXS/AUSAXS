#pragma once

namespace hist {
    class DistanceHistogram;
    class CompositeDistanceHistogram;
    class CompositeDistanceHistogramFF;
    class Histogram;
    class HistogramManager;

    /**
     * @brief A ScatteringProfile is just a (q, I(q)) histogram. 
     */    
    using ScatteringProfile = Histogram;
}