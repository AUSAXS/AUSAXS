// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>

namespace ausaxs::settings {
    // Rigid-body optimization is configured entirely by its script (or through the C++ API), so there are no global
    // settings here: each choice is an argument to the matching factory. Only the vocabulary remains, which the
    // factories and the script parsers share.
    struct EXPORT rigidbody {
        enum class TransformationStrategyChoice {
            RigidTransform,     // Transform all bodies connected to one side of the constraint. 
            SingleTransform,    // Transform only the body directly connected to one side of the constraint.
            ForceTransform      // Rotations and translations are applied as forces, resulting in more natural conformations. 
        };

        enum class BodySelectStrategyChoice {
            RandomBodySelect,           // Select a random body, then a random constraint within that body. 
            RandomConstraintSelect,     // Select a random constraint. 
            SequentialBodySelect,       // Select the first constraint, then the second, etc.
            SequentialConstraintSelect, // Select the first body, then the second, etc.
            ManualSelect                // Select a body and a constraint manually.
        };

        enum class ParameterGenerationStrategyChoice {
            Simple,             // Generate translation and rotation parameters.
            RotationsOnly,      // Only generate rotation parameters.
            TranslationsOnly,   // Only generate translation parameters.
            SymmetryOnly        // Only generate symmetry parameters.
        };

        enum class ParameterMaskStrategyChoice {
            All,                // All parameter components are active every step (default).
            Real,               // Only the real body transform (translation + rotation) is active.
            Symmetry,           // Only symmetry parameters (offset, axis, angle) are active.
            SymmetryTranslation,// Only the symmetry offset translation is active.
            SymmetryAxis,       // Only the symmetry frame orientation is active.
            Sequential,         // Alternates between real-only and symmetry-only steps.
            SequentialSymmetry, // Alternates between all-symmetry and symmetry-only steps, isolating the effect of the symmetry parameters.
            SequentialReal,     // Alternates between all-real and real-only steps, isolating the effect of the real body transform parameters.
            Random              // Randomly picks real-only or symmetry-only each step.
        };

        enum class DecayStrategyChoice {
            None,
            Linear,
            Exponential
        };

        enum class ConstraintGenerationStrategyChoice {
            None,       // Do not generate constraints. Only those supplied by the user will be used.
            Backbone    // Generate a bond constraint between every pair of backbone-adjacent bodies.
        };

        enum class ControllerChoice {
            Classic,    // Classic controller essentially equivalent to a gradient descent. 
            Metropolis, // Metropolis controller relying on Bayesian statistics to find probable conformations.
        };
    };
}