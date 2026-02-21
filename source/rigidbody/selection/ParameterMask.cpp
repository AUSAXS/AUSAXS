// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ParameterMask.h>

#include <cassert>

using namespace ausaxs::rigidbody::selection;

void ParameterMask::apply(parameter::BodyTransformParametersRelative& params) const {
    assert(
        (real_translation || real_rotation || sym_translation || sym_axis)
        && "ParameterMask::apply: all mask fields are false — nothing will be optimized."
    );

    if (!real_translation) {params.translation = std::nullopt;}
    if (!real_rotation)    {params.rotation    = std::nullopt;}

    if (params.symmetry_pars.has_value()) {
        if (!sym_translation && !sym_axis) {params.symmetry_pars = std::nullopt;} 
        else {
            for (auto& sym : params.symmetry_pars.value()) {
                if (!sym_translation) {for(auto& t : sym->span_translation()) {t = 0;}}
                if (!sym_axis) {for(auto& r : sym->span_rotation()) {r = 0;}}
            }
        }
    }
}
