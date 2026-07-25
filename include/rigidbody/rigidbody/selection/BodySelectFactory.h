// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace ausaxs::rigidbody {
    namespace factory {
        std::unique_ptr<selection::BodySelectStrategy> create_selection_strategy(observer_ptr<const Rigidbody> molecule);

        std::unique_ptr<selection::BodySelectStrategy> create_selection_strategy(
            observer_ptr<const Rigidbody> molecule, settings::rigidbody::BodySelectStrategyChoice choice
        );

        std::unique_ptr<selection::BodySelectStrategy> create_selection_strategy(
            observer_ptr<const Rigidbody> molecule,
            settings::rigidbody::BodySelectStrategyChoice body_choice,
            settings::rigidbody::ParameterMaskStrategyChoice mask_choice
        );

        /**
         * @brief Create a ManualSelect strategy that always selects the given body.
         *
         * @param ibody The index of the body to select on every call.
         */
        std::unique_ptr<selection::BodySelectStrategy> create_manual_selection_strategy(
            observer_ptr<const Rigidbody> molecule, unsigned int ibody
        );

        std::unique_ptr<selection::BodySelectStrategy> create_manual_selection_strategy(
            observer_ptr<const Rigidbody> molecule, unsigned int ibody, settings::rigidbody::ParameterMaskStrategyChoice mask_choice
        );

        /**
         * @brief Create a ManualSelect strategy that always selects a single declared symmetry of the given body.
         *
         * Only that symmetry's parameters are optimized; the host body's rigid pose and its other symmetries stay frozen.
         * Unlike the overloads above, this deliberately ignores settings::rigidbody::parameter_mask_strategy — targeting one
         * symmetry is an unconditional intent, not a default to be overridden.
         *
         * @param ibody The index of the body hosting the symmetry.
         * @param isymmetry The index of the symmetry within that body.
         */
        std::unique_ptr<selection::BodySelectStrategy> create_manual_symmetry_selection_strategy(
            observer_ptr<const Rigidbody> molecule, unsigned int ibody, unsigned int isymmetry
        );
    }
}