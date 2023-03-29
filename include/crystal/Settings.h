#pragma once

namespace setting {
    struct crystal {
        enum class MillerGenerationChoice {All, Fibonacci, Reduced};
        inline static MillerGenerationChoice mgc = MillerGenerationChoice::All; // The choice of Miller index generation algorithm.
        inline static unsigned int h = 100;                                     // The maximum Miller index along the x direction.
        inline static unsigned int k = 100;                                     // The maximum Miller index along the y direction.
        inline static unsigned int l = 100;                                     // The maximum Miller index along the z direction.

        inline static double max_q = 1e6;
    };
}