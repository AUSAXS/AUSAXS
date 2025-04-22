#pragma once

#include <md/programs/all.h>
#include <md/programs/mdrun/Execution.h>
#include <md/simulate/GMXOptions.h>
#include <md/utility/files/MDPCreator.h>
#include <utility/StringUtils.h>

namespace ausaxs::md {
    struct SimulateMoleculeOptions {
        SystemSettings system;
        std::string jobname;
        PDBFile pdbfile;
        MDPFile mdp;
        RunLocation setup_runner = RunLocation::local;
        RunLocation main_runner = RunLocation::local;
        std::string jobscript;
    };

    struct SimulateMoleculeOutput {
        std::shared_ptr<shell::Jobscript<MDRunResult>> job;        
        TOPFile top;    // topology file
        GROFile gro;    // solv_ion structure
    };

    SimulateMoleculeOutput simulate_molecule(SimulateMoleculeOptions&& options);
}