#pragma once

#include <string>

namespace ausaxs::settings::md {
    extern std::string gmx_path;                    // The path to the GROMACS installation.
    extern std::string buffer_path;                 // The path to the buffer directory.
    extern std::string water_model;                 // The water model to use.
    extern std::string force_field;                 // The force field to use.
    extern std::string box_type;                    // The geometric object to use as an envelope.
    extern std::string cation;                      // The cation atom to use.
    extern std::string anion;                       // The anion atom to use.
    extern std::string minimization_sim_location;   // The location of the minimization simulation.
    extern std::string thermalization_sim_location; // The location of the thermalization simulation.
    extern std::string production_sim_location;     // The location of the production simulation.
}