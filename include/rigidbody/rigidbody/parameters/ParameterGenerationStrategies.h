// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <utility/Random.h>

namespace ausaxs::rigidbody::parameter {
    template<bool TRANSLATE, bool ROTATE, bool SYMMETRY>
    class LimitedParameterGenerator : public ParameterGenerationStrategy {
        public: 
            using ParameterGenerationStrategy::ParameterGenerationStrategy;
            ~LimitedParameterGenerator() override = default;

            BodyTransformParameters next(int ibody) override;
    };
}

template<bool TRANSLATE, bool ROTATE, bool SYMMETRY>
ausaxs::rigidbody::parameter::BodyTransformParameters ausaxs::rigidbody::parameter::LimitedParameterGenerator<TRANSLATE, ROTATE, SYMMETRY>::next(int ibody) {
    double scaling = decay_strategy->next();

    assert(ibody < static_cast<int>(rigidbody->conformation->configuration.parameters.size()) && "ibody out of bounds");
    auto& current_pars = rigidbody->conformation->configuration.parameters[ibody];

    Vector3<double> t = current_pars.translation;
    if constexpr (TRANSLATE) {
        t.x() += translation_dist(random::generator())*scaling;
        t.y() += translation_dist(random::generator())*scaling;
        t.z() += translation_dist(random::generator())*scaling;
    }

    Vector3<double> r = current_pars.rotation;
    if constexpr (ROTATE) {
        r.x() += rotation_dist(random::generator())*scaling;
        r.y() += rotation_dist(random::generator())*scaling;
        r.z() += rotation_dist(random::generator())*scaling;
    }

    std::vector<symmetry::Symmetry> symmetry_pars = current_pars.symmetry_pars;
    if constexpr (SYMMETRY) {
        auto symmetries = static_cast<const ausaxs::symmetry::OptimizableSymmetryStorage*>(rigidbody->molecule.get_body(ibody).symmetry().get_obj());
        for (size_t i = 0; i < symmetries->symmetries.size(); ++i) {
            if (symmetries->optimize_translate) {
                symmetry_pars[i].initial_relation.translation.x() += symmetry_dist(random::generator())*scaling;
                symmetry_pars[i].initial_relation.translation.y() += symmetry_dist(random::generator())*scaling;
                symmetry_pars[i].initial_relation.translation.z() += symmetry_dist(random::generator())*scaling;
            }

            if (symmetries->optimize_rotate) {
                symmetry_pars[i].initial_relation.orientation.x() += symmetry_dist(random::generator())*scaling;
                symmetry_pars[i].initial_relation.orientation.y() += symmetry_dist(random::generator())*scaling;
                symmetry_pars[i].initial_relation.orientation.z() += symmetry_dist(random::generator())*scaling;
            }
        }
    }

    return {t, r, std::move(symmetry_pars)};
}

namespace ausaxs::rigidbody::parameter {
    using AllParameters    = LimitedParameterGenerator<true,  true,  true>;
    using TranslationsOnly = LimitedParameterGenerator<true,  false, false>;
    using RotationsOnly    = LimitedParameterGenerator<false, true,  false>;
    using SymmetryOnly     = LimitedParameterGenerator<false, false, true>;
}