#pragma once

#include <md/simulate/GMXOptions.h>
#include <md/programs/mdrun/MDExecutor.h>
#include <md/utility/files/all.h>

namespace ausaxs::md {
    struct SimulateBufferOptions {
        SystemSettings system;
        GROFile refgro;
        MDPFile mdp;
        std::unique_ptr<executor::type> minimize_runner;
        std::unique_ptr<executor::type> equilibrate_runner;
        std::unique_ptr<executor::type> production_runner;
    };

    struct SimulateBufferOutput {
        std::shared_ptr<Executor<MDRunResult>> job;
        TOPFile top;    // topology file
        GROFile gro;    // production structure
    };

    SimulateBufferOutput simulate_buffer(SimulateBufferOptions&& options);
}