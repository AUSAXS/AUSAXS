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
            auto& symmetries = params.symmetry_pars.value();
            assert(
                (!target_symmetry.has_value() || *target_symmetry < symmetries.size())
                && "ParameterMask::apply: target_symmetry is out of range for the generated parameters."
            );

            for (unsigned int i = 0; i < symmetries.size(); ++i) {
                // a targeted mask freezes every symmetry but the targeted one; an untargeted one treats them all alike
                bool active = !target_symmetry.has_value() || i == *target_symmetry;
                bool translate = sym_translation && active;
                bool rotate = sym_axis && active;

                // recurse into composite sub-symmetries; they have no contiguous span of their own
                symmetry::for_each_leaf(*symmetries[i], [&](symmetry::ISymmetry& leaf) {
                    if (!translate) {for (auto& t : leaf.span_translation()) {t = 0;}}
                    if (!rotate) {for (auto& r : leaf.span_rotation()) {r = 0;}}
                });
            }
        }
    }
}
