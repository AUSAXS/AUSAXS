// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/decay/ExponentialDecay.h>

#include <cmath>

using namespace ausaxs::rigidbody::parameter::decay;

ExponentialDecay::ExponentialDecay(unsigned int max_iterations) {
    set_characteristic_time(max_iterations/2);
}

ExponentialDecay::~ExponentialDecay() = default;

double ExponentialDecay::next() {
    return std::exp(-decay_rate*draws++);
}

void ExponentialDecay::set_characteristic_time(unsigned int iterations) {
    decay_rate = 1.0/iterations;
}