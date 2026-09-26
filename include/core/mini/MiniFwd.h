// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::mini {
    struct Parameter;
    struct FittedParameter;
    class Landscape;

    enum class algorithm { //NOLINT
        GOLDEN,
        MINIMUM_EXPLORER,
        SCAN,
        LIMITED_SCAN,
        LEVENBERG_MARQUARDT,
        #if defined(DLIB_AVAILABLE)
            DLIB_GLOBAL,
            BFGS,
        #endif
        DEFAULT=LEVENBERG_MARQUARDT
    };
}