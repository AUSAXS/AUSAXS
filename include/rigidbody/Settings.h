#pragma once

#include <vector>
#include <string>

// namespace setting {
//     struct rigidbody {
//         enum class TransformationStrategyChoice {
//             RigidTransform,     // Transform all bodies connected to one side of the constraint. 
//             SingleTransform,    // Transform only the body directly connected to one side of the constraint.
//             ForceTransform      // Rotations and translations are applied as forces, resulting in more natural conformations. 
//         };
//         enum class ParameterGenerationStrategyChoice {
//             Simple,         // Generate translation and rotation parameters. Their amplitudes decays linearly with the iteration number.
//             RotationsOnly   // Only generate rotation parameters. The amplitudes decays linearly with the iteration number.
//         };
//         enum class BodySelectStrategyChoice {
//             RandomSelect,           // Select a random body, then a random constraint within that body. 
//             RandomConstraintSelect, // Select a random constraint. 
//             SequentialSelect        // Select the first constraint, then the second, etc.
//         };

//         inline static TransformationStrategyChoice tsc = TransformationStrategyChoice::RigidTransform;
//         inline static ParameterGenerationStrategyChoice pgsc = ParameterGenerationStrategyChoice::Simple;
//         inline static BodySelectStrategyChoice bssc = BodySelectStrategyChoice::RandomSelect;

//         inline static unsigned int iterations = 1000;   // The number of iterations to run the rigid body optimization for.
//         inline static double bond_distance = 3;         // The maximum distance in Ångström between two atoms that allows for a constraint. 

//         struct detail {
//             inline static std::vector<int> constraints; // The residue ids to place a constraint at.
//             inline static std::string calibration_file; // The file to read constraints from.
//         };
//     };
// }