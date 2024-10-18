#pragma once

namespace mini {
    struct Parameter;
    struct FittedParameter;
    class Landscape;

    enum class algorithm {
        GOLDEN,
        MINIMUM_EXPLORER,
        SCAN,
        LIMITED_SCAN,
        #if defined(DLIB_AVAILABLE)
            DLIB_GLOBAL,
            BFGS,
            DEFAULT=BFGS
        #else
            DEFAULT=GOLDEN
        #endif
    };
}