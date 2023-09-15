#pragma once

namespace settings {
    namespace axes {
        extern unsigned int max_distance;  // The maximum distance in the p(r) function in Ångström. Should always be much larger than the actual maximum distance.
        extern double distance_bin_width;  // The width of each bin for the scattering plots.
        extern unsigned int bins;          // The number of bins for the scattering plots.
        extern double qmin;                // Lower limit on the used q-values
        extern double qmax;                // Upper limit on the used q-values
        extern unsigned int skip;          // The number of points to skip from the top of the scattering curve.
    }   
}

namespace settings::hist {
    enum class HistogramManagerChoice {
        HistogramManager,           // A simple manager that recalculates the entire histogram every time.
        HistogramManagerMT,         // A multithreaded implementation of the simple manager.
        HistogramManagerMTFF,       // A multithreaded implementation of the simple manager that uses precalculated form factor products.
        PartialHistogramManager,    // A smart manager that only recalculates the parts of the histogram that are needed.
        PartialHistogramManagerMT,  // A multithreaded implementation of the partial manager.
        PartialHistogramManagerMTFF // A multithreaded implementation of the partial manager that uses precalculated form factor products.
    };
    extern HistogramManagerChoice histogram_manager;
}