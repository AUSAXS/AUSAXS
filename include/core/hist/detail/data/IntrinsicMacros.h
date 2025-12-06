// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// AVX implies SSE4.1 and SSE2, but the MSVC compiler doesn't seem to define the latter two
#if defined __AVX__
    #if !defined __SSE2__
        #define __SSE2__
    #endif
    #if !defined __SSE4_1__
        #define __SSE4_1__
    #endif
#endif