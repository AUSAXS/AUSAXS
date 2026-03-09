// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// When AUSAXS_FORCE_SCALAR is defined, suppress all SIMD dispatch so the compiler
// auto-vectorizes the scalar fallback paths instead of using hand-written intrinsics.
// Used in per-SIMD benchmark variants to obtain a fair comparison baseline.
#if defined AUSAXS_FORCE_SCALAR
    #undef __AVX512F__
    #undef __AVX__
    #undef __SSE4_1__
    #undef __SSE2__
#else

// AVX implies SSE4.1 and SSE2, but the MSVC compiler doesn't seem to define the latter two
#if defined __AVX__
    #if !defined __SSE2__
        #define __SSE2__
    #endif
    #if !defined __SSE4_1__
        #define __SSE4_1__
    #endif
#endif

#endif // AUSAXS_FORCE_SCALAR