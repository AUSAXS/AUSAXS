#pragma once

namespace setting {
    struct crystal {
        enum class MillerGenerationChoice {All, Fibonacci, Reduced};
        inline static MillerGenerationChoice mgc = MillerGenerationChoice::All; // The choice of Miller index generation algorithm.
        inline static unsigned int h = 100;                                     // The maximum Miller index along the x direction.
        inline static unsigned int k = 100;                                     // The maximum Miller index along the y direction.
        inline static unsigned int l = 100;                                     // The maximum Miller index along the z direction.

        inline static double max_q = 1e6;                                       // The maximum length of the Miller indices. 
        inline static double grid_expansion = 3;                                // The factor by which the grid is expanded when loading a pdb structure. 

        struct reduced {
            inline static double basis_q = 3;                                   // The maximum q value for which the basis is generated.
        };
    };
}