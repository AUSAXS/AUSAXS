// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::detail {
    /**
     * @brief Diagnostic printers for the assert lambdas in the math headers.
     */
    bool report_index_1d(const char* context, long long i, long long bound);
    bool report_index_2d(const char* context, long long i, long long j, long long n, long long m);
    bool report_index_3d(const char* context, long long i, long long j, long long k, long long n, long long m, long long l);

    /**
     * @brief Report a mismatch between a single observed and expected extent.
     */
    bool report_size_mismatch(const char* context, long long got, long long expected);

    /**
     * @brief Report a mismatch between an observed and expected pair of extents.
     */
    bool report_dim_mismatch(const char* context, long long got_n, long long got_m, long long expected_n, long long expected_m);
}