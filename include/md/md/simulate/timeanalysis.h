#pragma once

#include <md/simulate/GMXOptions.h>

namespace md {
    /**
     * @brief Generate a series of SAXS curves with different durations. 
     * 
     * @param options The run settings. 
     * @param timestep The additional duration for each run in picoseconds.
     */
    std::vector<SAXSOutput> timeanalysis(SAXSOptions& options, double timestep);
}