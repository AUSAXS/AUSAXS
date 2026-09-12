// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <math/detail/Diagnostics.h>

#include <iostream>

using namespace ausaxs;

bool detail::report_index_1d(const char* context, long long i, long long bound) {
    std::cout << context << ": Index out of bounds (" << i << " should be less than " << bound << ")" << std::endl;
    return false;
}

bool detail::report_index_2d(const char* context, long long i, long long j, long long n, long long m) {
    std::cout << context << ": Index out of bounds (" << i << ", " << j << ") should be less than (" << n << ", " << m << ")" << std::endl;
    return false;
}

bool detail::report_index_3d(const char* context, long long i, long long j, long long k, long long n, long long m, long long l) {
    std::cout << context << ": Index out of bounds (" << i << ", " << j << ", " << k << ") should be less than (" << n << ", " << m << ", " << l << ")" << std::endl;
    return false;
}

bool detail::report_size_mismatch(const char* context, long long got, long long expected) {
    std::cout << context << ": Dimensions do not match (got: " << got << ", expected: " << expected << ")." << std::endl;
    return false;
}

bool detail::report_dim_mismatch(const char* context, long long got_n, long long got_m, long long expected_n, long long expected_m) {
    std::cout << context << ": Dimensions do not match (got: [" << got_n << ", " << got_m << "], expected: [" << expected_n << ", " << expected_m << "])." << std::endl;
    return false;
}
