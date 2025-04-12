#pragma once

#include <md/simulate/saxs.h>

namespace ausaxs::md {
    /**
     * @brief Generate a series of SAXS curves with different durations. 
     * 
     * @param options The run settings. 
     * @param timestep The additional duration for each run in picoseconds.
     */
    std::vector<SAXSOutput> timeanalysis(SimulateSAXSOptions&& options, double timestep);
}