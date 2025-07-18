// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::mini {
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
            DEFAULT=DLIB_GLOBAL
        #else
            DEFAULT=GOLDEN
        #endif
    };
}