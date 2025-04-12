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

    SimulateBufferOutput simulate_buffer(SimulateBufferOptions&& options);
}