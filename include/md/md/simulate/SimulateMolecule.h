#pragma once

#include <md/programs/all.h>
#include <md/programs/mdrun/MDExecutor.h>
#include <md/simulate/GMXOptions.h>
#include <md/utility/files/MDPCreator.h>
#include <utility/StringUtils.h>

namespace ausaxs::md {
    struct SimulateMoleculeOptions {
        SystemSettings system;
        PDBFile pdbfile;
        MDPFile mdp_equilibration;
        MDPFile mdp_production;
        std::unique_ptr<executor::type> minimize_runner;
        std::unique_ptr<executor::type> equilibrate_runner;
        std::unique_ptr<executor::type> production_runner;
    };

    struct SimulateMoleculeOutput {
        std::shared_ptr<Executor<MDRunResult>> job;
        TOPFile top;    // topology file
        GROFile gro;    // solv_ion structure
    };

    SimulateMoleculeOutput simulate_molecule(SimulateMoleculeOptions&& options);
}