// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/SettingRef.h>
#include <settings/SettingsIORegistry.h>
#include <settings/ExportMacro.h>

#include <vector>
#include <string>

namespace ausaxs::settings {
    struct EXPORT rigidbody {
        static unsigned int iterations;   // The number of iterations to run the rigid body optimization for.
        static double bond_distance;      // The maximum distance in Ångström between two atoms that allows for a constraint.

        struct detail {
            static std::vector<int> constraints; // The residue ids to place a constraint at.
            static std::string calibration_file; // The file to read constraints from.
        };

        enum class TransformationStrategyChoice {
            RigidTransform,     // Transform all bodies connected to one side of the constraint. 
            SingleTransform,    // Transform only the body directly connected to one side of the constraint.
            ForceTransform      // Rotations and translations are applied as forces, resulting in more natural conformations. 
        };
        static TransformationStrategyChoice transform_strategy;

        enum class BodySelectStrategyChoice {
            RandomBodySelect,           // Select a random body, then a random constraint within that body. 
            RandomConstraintSelect,     // Select a random constraint. 
            SequentialBodySelect,       // Select the first constraint, then the second, etc.
            SequentialConstraintSelect, // Select the first body, then the second, etc.
            ManualSelect                // Select a body and a constraint manually.
        };
        static BodySelectStrategyChoice body_select_strategy;

        enum class ParameterGenerationStrategyChoice {
            Simple,             // Generate translation and rotation parameters.
            RotationsOnly,      // Only generate rotation parameters.
            TranslationsOnly,   // Only generate translation parameters.
            SymmetryOnly        // Only generate symmetry parameters.
        };
        static ParameterGenerationStrategyChoice parameter_generation_strategy;

        enum class DecayStrategyChoice {
            None,
            Linear,
            Exponential
        };
        static DecayStrategyChoice decay_strategy;

        enum class ConstraintGenerationStrategyChoice {
            None,       // Do not generate constraints. Only those supplied by the user will be used.
            Linear,     // Generate a linear chain of constraints between bodies.
            Volumetric  // Generate constraints between bodies based on proximity. 
        };
        static ConstraintGenerationStrategyChoice constraint_generation_strategy;
    };
}