// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <math/MovingAverager.h>

#include <math/Exceptions.h>

using namespace ausaxs;

void MovingAverage::validate_input(int N, int window_size) {
    if (N < window_size) {
        throw math::except::invalid_argument("MovingAverager::average: Window size is larger than data size.");
    }

    if (window_size % 2 == 0) {
        throw math::except::invalid_argument("Error in MovingAverager::validate_input: Window_size must be odd");
    }
}
