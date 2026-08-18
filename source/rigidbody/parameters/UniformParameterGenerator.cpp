// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <data/symmetry/CompositeSymmetry.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

BodyTransformParametersRelative UniformParameterGenerator::next(int ibody) {
    double scaling = decay_strategy->next();
    assert(ibody < static_cast<int>(rigidbody->conformation->absolute_parameters.parameters.size()) && "ibody out of bounds");
    BodyTransformParametersRelative params;

    if (amplitudes.translation != 0) {
        params.translation = {
            draw(amplitudes.translation)*scaling,
            draw(amplitudes.translation)*scaling,
            draw(amplitudes.translation)*scaling
        };
    }

    if (amplitudes.rotation != 0) {
        params.rotation = {
            draw(amplitudes.rotation)*scaling,
            draw(amplitudes.rotation)*scaling,
            draw(amplitudes.rotation)*scaling
        };
    }

    if (amplitudes.symmetry_translation != 0 || amplitudes.symmetry_rotation != 0) {
        // ReferenceSymmetryViews among these are inert: their spans are empty, so the loops below simply do not touch them.
        const auto& symmetries = rigidbody->molecule.get_body(ibody).symmetry().get_obj()->symmetries;
        params.symmetry_pars.emplace();

        // assign a random delta to every leaf sub-symmetry. CompositeSymmetry has no parameter spans of its own (its two 
        // sets cannot be expressed as one contiguous span), so we reach the leaves via for_each_leaf and perturb each of them.
        // a disabled component is explicitly zeroed rather than left alone, since the delta starts out as a clone of the symmetry itself.
        bool translate = amplitudes.symmetry_translation != 0;
        bool rotate = amplitudes.symmetry_rotation != 0;
        auto randomize = [&](ausaxs::symmetry::ISymmetry& leaf) {
            for (auto& t : leaf.span_translation()) {t = translate ? draw(amplitudes.symmetry_translation)*scaling : 0;}
            for (auto& r : leaf.span_rotation())    {r = rotate    ? draw(amplitudes.symmetry_rotation)*scaling    : 0;}
        };

        for (const auto& symmetry : symmetries) {
            auto delta = symmetry->clone();
            ausaxs::symmetry::for_each_leaf(*delta, randomize);
            params.symmetry_pars->emplace_back(std::move(delta));
        }
    }

    return params;
}
