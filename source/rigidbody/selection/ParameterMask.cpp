// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ParameterMask.h>
#include <data/symmetry/CompositeSymmetry.h>

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
                // recurse into composite sub-symmetries; they have no contiguous span of their own
                symmetry::for_each_leaf(*sym, [&](symmetry::ISymmetry& leaf) {
                    if (!sym_translation) {for (auto& t : leaf.span_translation()) {t = 0;}}
                    if (!sym_axis) {for (auto& r : leaf.span_rotation()) {r = 0;}}
                });
            }
        }
    }
}
