// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ParameterMask.h>

#include <cassert>

using namespace ausaxs::rigidbody::selection;

void ParameterMask::apply(parameter::BodyTransformParametersRelative& params) const {
    assert(
        (real_translation || real_rotation || sym_translation || sym_axis || sym_angle)
        && "ParameterMask::apply: all mask fields are false — nothing will be optimized."
    );

    if (!real_translation) {params.translation = std::nullopt;}
    if (!real_rotation)    {params.rotation    = std::nullopt;}

    if (params.symmetry_pars.has_value()) {
        if (!sym_translation && !sym_axis && !sym_angle) {
            params.symmetry_pars = std::nullopt;
        } else {
            for (auto& sym : params.symmetry_pars.value()) {
                if (!sym_translation) {sym->initial_relation.translation = {0, 0, 0};}
                if (!sym_axis) {sym->repeat_relation.axis = {0, 0, 0};}
                if (!sym_angle) {sym->repeat_relation.angle = 0;}
            }
        }
    }
}
