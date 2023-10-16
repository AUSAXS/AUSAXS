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
        HistogramManager,           // A simple manager that recalculates the entire histogram every time.
        HistogramManagerMT,         // A multithreaded implementation of the simple manager.
        HistogramManagerMTFF,       // A multithreaded implementation of the simple manager that uses precalculated form factor products.
        PartialHistogramManager,    // A smart manager that only recalculates the parts of the histogram that are needed.
        PartialHistogramManagerMT,  // A multithreaded implementation of the partial manager.
        PartialHistogramManagerMTFF,// A multithreaded implementation of the partial manager that uses precalculated form factor products.
        DebugManager
    };
    extern HistogramManagerChoice histogram_manager;
}