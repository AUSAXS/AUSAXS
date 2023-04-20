#pragma once

#include <utility/SmartOption.h>

namespace settings {
    namespace fit {
        extern settings::detail::SmartOption<bool> verbose;                 // Decides if the fitting process will be verbose.
        extern settings::detail::SmartOption<unsigned int> N;               // Number of points sampled when discretizing a model scattering curve
        extern settings::detail::SmartOption<unsigned int> max_iterations;  // Maximum number of iterations in the fitting process
    }
}