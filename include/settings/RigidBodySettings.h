#pragma once

#include <settings/SettingRef.h>
#include <settings/SettingsIORegistry.h>

#include <vector>
#include <string>

namespace settings {
    namespace rigidbody {
        extern unsigned int iterations;   // The number of iterations to run the rigid body optimization for.
        extern double bond_distance;      // The maximum distance in Ångström between two atoms that allows for a constraint.

        namespace detail {
            extern std::vector<int> constraints; // The residue ids to place a constraint at.
            extern std::string calibration_file; // The file to read constraints from.
        }
    }
}

namespace settings::rigidbody {
    enum class TransformationStrategyChoice {
        RigidTransform,     // Transform all bodies connected to one side of the constraint. 
        SingleTransform,    // Transform only the body directly connected to one side of the constraint.
        ForceTransform      // Rotations and translations are applied as forces, resulting in more natural conformations. 
    };
    extern TransformationStrategyChoice transform_strategy;
}

namespace settings::rigidbody {
    enum class BodySelectStrategyChoice {
        RandomSelect,           // Select a random body, then a random constraint within that body. 
        RandomConstraintSelect, // Select a random constraint. 
        SequentialSelect        // Select the first constraint, then the second, etc.
    };
    extern BodySelectStrategyChoice body_select_strategy;
}

namespace settings::rigidbody {
    enum class ParameterGenerationStrategyChoice {
        Simple,             // Generate translation and rotation parameters.
        RotationsOnly,      // Only generate rotation parameters.
        TranslationsOnly    // Only generate translation parameters.
    };
    extern ParameterGenerationStrategyChoice parameter_generation_strategy;
}

namespace settings::rigidbody {
    enum class DecayStrategyChoice {
        Linear,
        Exponential
    };
    extern DecayStrategyChoice decay_strategy;
}

namespace settings::rigidbody {
    enum class ConstraintGenerationStrategyChoice {
        None,       // Do not generate constraints. Only those supplied by the user will be used.
        Linear,     // Generate a linear chain of constraints between bodies.
        Volumetric  // Generate constraints between bodies based on proximity. 
    };
    extern ConstraintGenerationStrategyChoice constraint_generation_strategy;
}