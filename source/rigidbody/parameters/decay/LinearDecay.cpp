// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/decay/LinearDecay.h>

using namespace ausaxs::rigidbody::parameter::decay;

LinearDecay::LinearDecay(unsigned int max_iterations) : DecayStrategy(max_iterations) {
    set_characteristic_time(max_iterations/2);
}

LinearDecay::~LinearDecay() = default;

double LinearDecay::next() {
    return 1.0 - decay_rate*draws++;
}

void LinearDecay::set_characteristic_time(unsigned int iterations) {
    decay_rate = 0.5/iterations;
}