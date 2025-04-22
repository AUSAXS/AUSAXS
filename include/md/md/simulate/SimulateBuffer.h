#pragma once

#include <md/simulate/GMXOptions.h>

namespace ausaxs::md {
    struct SimulateBufferOptions {
        SystemSettings system;
        std::string jobname;
        GROFile refgro;
        MDPFile mdp;
        RunLocation setup_runner = RunLocation::local;
        RunLocation main_runner = RunLocation::local;
        std::string jobscript;
    };

    struct SimulateBufferOutput {
        std::shared_ptr<shell::Jobscript<MDRunResult>> job;
        TOPFile top;    // topology file
        GROFile gro;    // production structure
    };

    SimulateBufferOutput simulate_buffer(SimulateBufferOptions&& options);
}