// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// AVX2 implies SSE4.1 and SSE2, but MSVC doesn't always define the latter two
#if defined __AVX2__
    #if !defined __SSE2__
        #define __SSE2__
    #endif
    #if !defined __SSE4_1__
        #define __SSE4_1__
    #endif
#endif

// Resolve the effective SIMD dispatch level into AUSAXS_USE_{SSE2,AVX2,AVX512}.
//
// Supported baseline targets: x86-64-v1/v2 (SSE2), x86-64-v3 (AVX2+FMA), x86-64-v4 (AVX-512).
// Benchmark builds define one of USE_SCALAR, USE_SSE2, USE_AVX2, USE_AVX512
// to set a ceiling on dispatch, acting as safeguards alongside compiler flags.
// On Apple, SIMD dispatch is suppressed entirely (Clang compatibility).
// Otherwise, auto-detect from compiler builtins.
#if defined __APPLE__ || defined USE_SCALAR
    // No SIMD dispatch
#elif defined USE_SSE2
    #if defined __SSE2__
        #define AUSAXS_USE_SSE2
    #endif
#elif defined USE_AVX2
    #if defined __SSE2__
        #define AUSAXS_USE_SSE2
    #endif
    #if defined __AVX2__
        #define AUSAXS_USE_AVX2
    #endif
#elif defined USE_AVX512
    #if defined __SSE2__
        #define AUSAXS_USE_SSE2
    #endif
    #if defined __AVX2__
        #define AUSAXS_USE_AVX2
    #endif
    #if defined __AVX512F__
        #define AUSAXS_USE_AVX512
    #endif
#else
    // Auto-detect from compiler builtins
    #if defined __SSE2__
        #define AUSAXS_USE_SSE2
    #endif
    #if defined __AVX2__
        #define AUSAXS_USE_AVX2
    #endif
    #if defined __AVX512F__
        #define AUSAXS_USE_AVX512
    #endif
#endif