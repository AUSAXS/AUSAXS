#pragma once

#include <md/simulate/GMXOptions.h>

namespace ausaxs::md {
    struct SimulateSAXSOptions {
        std::string jobname;
        PDBFile pdbfile;
        SimulateMoleculeOutput molecule;
        SimulateBufferOutput buffer;
        RunLocation runner = RunLocation::local;
        std::string jobscript;
    };

    SAXSOutput simulate_saxs(SimulateSAXSOptions&& options);
}