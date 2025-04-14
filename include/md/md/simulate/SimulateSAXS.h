#pragma once

#include <md/simulate/SimulateMolecule.h>
#include <md/simulate/SimulateBuffer.h>

#include <optional>

namespace ausaxs::md {
    struct SimulateSAXSOptions {
        std::string jobname;
        PDBFile pdbfile;
        SimulateMoleculeOutput molecule;
        SimulateBufferOutput buffer;
        RunLocation runner = RunLocation::local;
        std::string jobscript;
        std::string output; 
        std::optional<mdp::templates::saxs::mol> mol_mdp;
        std::optional<mdp::templates::saxs::solv> buf_mdp;
    };

    struct SimulateSAXSOutput {
        std::unique_ptr<shell::Jobscript<SAXSRunResult>> job;
    };

    SimulateSAXSOutput simulate_saxs(SimulateSAXSOptions&& options);
}