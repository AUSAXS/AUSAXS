#pragma once

#include <md/simulate/saxs.h>

namespace ausaxs::md {
    std::vector<SAXSOutput> frameanalysis(SimulateSAXSOptions&& options);
}