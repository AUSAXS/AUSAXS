// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterAmplitudes.h>
#include <rigidbody/Rigidbody.h>

#include <numbers>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

double rigidbody::parameter::default_rotation() {
    return std::numbers::pi/3;
}

double rigidbody::parameter::default_symmetry_translation(observer_ptr<const Rigidbody> rigidbody) {
    return 2*rigidbody->molecule.get_Rg(false);
}

ParameterAmplitudes rigidbody::parameter::default_amplitudes(observer_ptr<const Rigidbody> rigidbody) {
    return {
        .translation = default_translation(),
        .rotation = default_rotation(),
        .symmetry_translation = default_symmetry_translation(rigidbody),
        .symmetry_rotation = default_symmetry_rotation()
    };
}
