// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterGenerationStrategies.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <utility/Random.h>

#include <functional>

template<bool TRANSLATE, bool ROTATE, bool SYMMETRY>
ausaxs::rigidbody::parameter::BodyTransformParametersRelative ausaxs::rigidbody::parameter::LimitedParameterGenerator<TRANSLATE, ROTATE, SYMMETRY>::next(int ibody) {
    double scaling = decay_strategy->next();
    assert(ibody < static_cast<int>(rigidbody->conformation->absolute_parameters.parameters.size()) && "ibody out of bounds");
    BodyTransformParametersRelative params;

    if constexpr (TRANSLATE) {
        params.translation = {
            translation_dist(random::generator())*scaling,
            translation_dist(random::generator())*scaling,
            translation_dist(random::generator())*scaling
        };
    }

    if constexpr (ROTATE) {
        params.rotation = {
            rotation_dist(random::generator())*scaling,
            rotation_dist(random::generator())*scaling,
            rotation_dist(random::generator())*scaling
        };
    }

    if constexpr (SYMMETRY) {
        auto symmetries = static_cast<const ausaxs::symmetry::OptimizableSymmetryStorage*>(rigidbody->molecule.get_body(ibody).symmetry().get_obj());
        params.symmetry_pars.emplace();

        // assign a random delta to every leaf sub-symmetry. CompositeSymmetry has no parameter
        // spans of its own (its two sets cannot be expressed as one contiguous span), so we
        // descend into its inner/outer parts and perturb each of them.
        std::function<void(ausaxs::symmetry::ISymmetry&)> randomize = [&](ausaxs::symmetry::ISymmetry& sym) {
            if (auto* comp = dynamic_cast<ausaxs::symmetry::CompositeSymmetry*>(&sym)) {
                randomize(*comp->inner);
                randomize(*comp->outer);
                return;
            }
            if (symmetries->optimize_translate) {
                for (auto& t : sym.span_translation()) {t = translation_symmetry_dist(random::generator())*scaling;}
            } else {for (auto& t : sym.span_translation()) {t = 0;}}

            if (symmetries->optimize_rot_axis) {
                for (auto& r : sym.span_rotation()) {r = rotation_symmetry_dist(random::generator())*scaling;}
            } else {for (auto& r : sym.span_rotation()) {r = 0;}}
        };

        for (int i = 0; i < static_cast<int>(symmetries->symmetries.size()); ++i) {
            auto delta = symmetries->symmetries[i]->clone();
            randomize(*delta);
            params.symmetry_pars->emplace_back(std::move(delta));
        }
    }

    return params;
}

template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<true,  true,  true>;
template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<true,  false, false>;
template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<false, true,  false>;
template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<false, false, true>;