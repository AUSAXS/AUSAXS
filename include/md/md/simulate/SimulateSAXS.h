#pragma once

#include <md/simulate/SimulateMolecule.h>
#include <md/simulate/SimulateBuffer.h>

namespace ausaxs::md {
    struct SimulateSAXSOptions {
        std::string jobname;
        PDBFile pdbfile;
        SimulateMoleculeOutput molecule;
        SimulateBufferOutput buffer;
        RunLocation runner = RunLocation::local;
        std::string jobscript;
    };

    struct SimulateSAXSOutput {
        std::unique_ptr<shell::Jobscript<SAXSRunResult>> job;
    };

    SimulateSAXSOutput simulate_saxs(SimulateSAXSOptions&& options);
}