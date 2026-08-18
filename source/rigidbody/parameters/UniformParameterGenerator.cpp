// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <math/Vector3.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

BodyTransformParametersRelative UniformParameterGenerator::next(int ibody) {
    double scaling = decay_strategy->next();
    assert(ibody < static_cast<int>(rigidbody->conformation->absolute_parameters.parameters.size()) && "ibody out of bounds");
    BodyTransformParametersRelative params;

    // every amplitude bounds the magnitude of the whole step, which is then pointed in a uniformly drawn direction. Drawing each
    // component independently instead would fill a cube, whose corners reach sqrt(3) further than its faces: both a larger step
    // than the amplitude names, and one biased towards the diagonals.
    if (amplitudes.translation != 0) {
        params.translation = draw_direction()*(draw(amplitudes.translation)*scaling);
    }

    if (amplitudes.rotation != 0) {
        params.rotation = draw_direction()*(draw(amplitudes.rotation)*scaling);
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
            // the rotation amplitude is an angle in radians whatever the symmetry stores internally, so the symmetry is asked to
            // express it in its own parameters. This is read off the leaf before anything is written to it, since the conversion
            // may depend on the leaf's current rotation parameters.
            auto rotation = rotate
                ? leaf.rotation_from_angle(draw(amplitudes.symmetry_rotation)*scaling, draw_direction())
                : Vector3<double>{0, 0, 0};

            auto t = leaf.span_translation();
            if (translate) {draw_isotropic(t, draw(amplitudes.symmetry_translation)*scaling);}
            else {for (auto& v : t) {v = 0;}}

            auto r = leaf.span_rotation();
            assert((r.empty() || r.size() == 3) && "UniformParameterGenerator::next: unexpected rotation parameter count.");
            for (std::size_t i = 0; i < r.size(); ++i) {r[i] = rotation[i];}
        };

        for (const auto& symmetry : symmetries) {
            auto delta = symmetry->clone();
            ausaxs::symmetry::for_each_leaf(*delta, randomize);
            params.symmetry_pars->emplace_back(std::move(delta));
        }
    }

    return params;
}
