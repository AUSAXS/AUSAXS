#pragma once

#include <md/simulate/SimulateMolecule.h>
#include <md/simulate/SimulateBuffer.h>

#include <optional>

namespace ausaxs::md {
    struct SimulateSAXSOptions {
        PDBFile pdbfile;
        SimulateMoleculeOutput molecule;
        SimulateBufferOutput buffer;
        std::string output; 
        std::unique_ptr<executor::type> runner;
        std::optional<mdp::templates::saxs::mol> mol_mdp;
        std::optional<mdp::templates::saxs::solv> buf_mdp;
    };

    struct SimulateSAXSOutput {
        std::unique_ptr<Executor<SAXSRunResult>> job;
    };

    SimulateSAXSOutput simulate_saxs(SimulateSAXSOptions&& options);
}