#pragma once

#include <utility/SmartOption.h>

#include <vector>
#include <string>

namespace setting {
    namespace rigidbody {
        extern settings::detail::SmartOption<unsigned int> iterations;   // The number of iterations to run the rigid body optimization for.
        extern settings::detail::SmartOption<double> bond_distance;      // The maximum distance in Ångström between two atoms that allows for a constraint.

        namespace detail {
            extern settings::detail::SmartOption<std::vector<int>> constraints; // The residue ids to place a constraint at.
            extern settings::detail::SmartOption<std::string> calibration_file; // The file to read constraints from.
        }
    }
}