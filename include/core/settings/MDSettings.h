#pragma once

#include <string>

namespace ausaxs::settings {
    struct md {
        static std::string gmx_path;                    // The path to the GROMACS installation.
        static std::string buffer_path;                 // The path to the buffer directory.
        static std::string water_model;                 // The water model to use.
        static std::string force_field;                 // The force field to use.
        static std::string box_type;                    // The geometric object to use as an envelope.
        static std::string cation;                      // The cation atom to use.
        static std::string anion;                       // The anion atom to use.
        static std::string minimization_sim_location;   // The location of the minimization simulation.
        static std::string thermalization_sim_location; // The location of the thermalization simulation.
        static std::string production_sim_location;     // The location of the production simulation.
    };
}