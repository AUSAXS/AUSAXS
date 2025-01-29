#pragma once

#include <settings/ExportMacro.h>

namespace ausaxs::settings {
    struct EXPORT axes {
        static unsigned int skip;   // The number of points to skip from the top of the scattering curve.
        static double qmin;         // Lower limit on the used q-values
        static double qmax;         // Upper limit on the used q-values
    };

    struct EXPORT hist {
        enum class HistogramManagerChoice {
            None,                                // No histogram manager. This is only used to avoid performing an expensive constructor when the manager is not needed.
            HistogramManager,                    // A simple manager that recalculates the entire histogram every time.
            HistogramManagerMT,                  // A multithreaded implementation of the simple manager.
            HistogramManagerMTFFAvg,             // A multithreaded implementation of the simple manager that uses precalculated form factor products and an average for the excluded volume.
            HistogramManagerMTFFExplicit,        // A multithreaded implementation of the simple manager that uses precalculated form factor products for both the protein and the excluded volume. 
            HistogramManagerMTFFGrid,            // A multithreaded implementation of the simple manager using a grid-based approach to evaluate the excluded volume. 
            HistogramManagerMTFFGridSurface,     // A multithreaded implementation of the simple manager using a grid-based approach to evaluate the excluded volume with a surface correction.
            HistogramManagerMTFFGridScalableExv, // A multithreaded implementation of the simple manager using a grid-based approach to evaluate a scalable excluded volume.
            HistogramSymmetryManagerMT,          // A multithreaded implementation of the partial symmetry manager.
            PartialHistogramManager,             // A smart manager that only recalculates the parts of the histogram that have been changed between each call. 
            PartialHistogramManagerMT,           // A multithreaded implementation of the partial manager.
            PartialHistogramManagerMTFFAvg,      // A multithreaded implementation of the partial manager that uses precalculated form factor products and an average for the excluded volume.
            PartialHistogramManagerMTFFExplicit, // A multithreaded implementation of the partial manager that uses precalculated form factor products for both the protein and the excluded volume. 
            PartialHistogramManagerMTFFGrid,     // A multithreaded implementation of the partial manager using a grid-based approach to evaluate the excluded volume.
            PartialHistogramSymmetryManagerMT,   // A multithreaded implementation of the partial symmetry manager.
            FoXSManager,                         // A manager that mimics the FoXS method to evaluate the scattering intensity.
            PepsiManager,                        // A manager that mimics the Pepsi method to evaluate the scattering intensity.
            CrysolManager,                       // A manager that mimics the Crysol method to evaluate the scattering intensity.
        };
        static bool weighted_bins;          // Whether to use weighted p(r) bins or not.
        static HistogramManagerChoice histogram_manager;
    };
}