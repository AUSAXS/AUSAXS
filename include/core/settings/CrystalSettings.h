#pragma once

namespace settings {
    namespace crystal {
        extern unsigned int h;        // The maximum Miller index along the x direction.
        extern unsigned int k;        // The maximum Miller index along the y direction.
        extern unsigned int l;        // The maximum Miller index along the z direction.

        extern double max_q;          // The maximum length of the Miller indices. 
        extern double grid_expansion; // The factor by which the grid is expanded when loading a pdb structure. 

        namespace reduced {
            extern double basis_q;    // The maximum q value for which the basis is generated.
        }

        namespace detail {
            extern bool use_checkpointing; // Whether to use checkpointing during the calculation. 
        }
    }
}

namespace settings::crystal {
    // The choice of Miller index generation algorithm. See the documentation for the individual classes for more information.
    enum class MillerGenerationChoice {
        All, // Generates all miller indices within the specified range. 
        Reduced, // Generates a subset of the miller indices within the specified range using a reduced basis.
        Fibonacci // Similar to Reduced, but uses a Fibonacci sphere to make the basis more uniform. Experimental and not recommended for general usage.
    };
    extern MillerGenerationChoice miller_generation_strategy;
}