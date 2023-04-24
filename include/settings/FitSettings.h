#pragma once

namespace settings {
    namespace fit {
        extern bool verbose;                 // Decides if the fitting process will be verbose.
        extern unsigned int N;               // Number of points sampled when discretizing a model scattering curve
        extern unsigned int max_iterations;  // Maximum number of iterations in the fitting process
    }
}