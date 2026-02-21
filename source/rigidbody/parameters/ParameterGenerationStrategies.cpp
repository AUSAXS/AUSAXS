// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/ParameterGenerationStrategies.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <utility/Random.h>

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
        params.symmetry_pars = std::vector<std::unique_ptr<symmetry::ISymmetry>>(symmetries->symmetries.size());
        for (int i = 0; i < static_cast<int>(params.symmetry_pars->size()); ++i) {
            auto& current_sym = params.symmetry_pars.value()[i];
            if (symmetries->optimize_translate) {
                current_sym->initial_relation.translation.x() = translation_symmetry_dist(random::generator())*scaling;
                current_sym->initial_relation.translation.y() = translation_symmetry_dist(random::generator())*scaling;
                current_sym->initial_relation.translation.z() = translation_symmetry_dist(random::generator())*scaling;
            }

            if (symmetries->optimize_rot_axis) {
                current_sym->repeat_relation.axis.x() = rotation_symmetry_dist(random::generator())*scaling;
                current_sym->repeat_relation.axis.y() = rotation_symmetry_dist(random::generator())*scaling;
                current_sym->repeat_relation.axis.z() = rotation_symmetry_dist(random::generator())*scaling;
            }

            if (symmetries->optimize_rot_angle) {
                current_sym->repeat_relation.angle = angle_symmetry_dist(random::generator())*scaling;
            }
        }
    }

    return params;
}

template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<true,  true,  true>;
template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<true,  false, false>;
template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<false, true,  false>;
template class ausaxs::rigidbody::parameter::LimitedParameterGenerator<false, false, true>;