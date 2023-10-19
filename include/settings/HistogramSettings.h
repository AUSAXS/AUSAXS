#pragma once

namespace settings {
    namespace axes {
        extern unsigned int skip;   // The number of points to skip from the top of the scattering curve.
        extern double qmin;         // Lower limit on the used q-values
        extern double qmax;         // Upper limit on the used q-values
    }   
}

namespace settings::hist {
    enum class HistogramManagerChoice {
        HistogramManager,                    // A simple manager that recalculates the entire histogram every time.
        HistogramManagerMT,                  // A multithreaded implementation of the simple manager.
        HistogramManagerMTFFAvg,             // A multithreaded implementation of the simple manager that uses precalculated form factor products and an average for the excluded volume.
        HistogramManagerMTFFExplicit,        // A multithreaded implementation of the simple manager that uses precalculated form factor products for both the protein and the excluded volume. 
        PartialHistogramManager,             // A smart manager that only recalculates the parts of the histogram that have been changed between each call. 
        PartialHistogramManagerMT,           // A multithreaded implementation of the partial manager.
        PartialHistogramManagerMTFFAvg,      // A multithreaded implementation of the partial manager that uses precalculated form factor products and an average for the excluded volume.
        PartialHistogramManagerMTFFExplicit, // A multithreaded implementation of the partial manager that uses precalculated form factor products for both the protein and the excluded volume. 
        DebugManager
    };
    extern HistogramManagerChoice histogram_manager;
}