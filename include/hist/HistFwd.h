#pragma once

namespace hist {
    class DebyeLookupTable;
    class DistanceHistogram;
    class CompositeDistanceHistogram;
    class CompositeDistanceHistogramFF;
    class Histogram;

    /**
     * @brief A ScatteringHistogram is just a (q, I(q)) histogram. 
     */    
    using ScatteringHistogram = Histogram;
}