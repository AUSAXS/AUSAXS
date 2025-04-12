#pragma once

#include <md/simulate/SimulateSAXS.h>

namespace ausaxs::md {
    std::vector<SimulateSAXSOutput> frameanalysis(SimulateSAXSOptions&& options);
}