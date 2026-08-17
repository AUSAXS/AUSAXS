// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::constraints {
    class DistanceConstraintBond : public DistanceConstraintAtom {
        public: 
            /**
             * @brief Create a backbone bond constraint between the two bodies.
             *
             * The representative C-alpha pair is selected automatically: only the terminal C-alphas of each body are eligible. 
             * Throws if no suitable pair is found or if the selected atoms are too far apart to represent a bond.
             */
            DistanceConstraintBond(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2);

            /**
             * @brief Restore a bond from stored values, skipping the candidate search entirely. See @ref restore_t.
             *
             * Notably this also skips the "too far apart" rejection: a refinement is free to stretch a bond well past
             * the 4 Å a *new* bond is allowed to span, and a state that was legal to write must stay legal to read.
             */
            DistanceConstraintBond(
                restore_t, observer_ptr<const data::Molecule> molecule,
                int ibody1, int iatom1, int ibody2, int iatom2,
                std::pair<int, int> isym1, std::pair<int, int> isym2, double d_target
            );

            virtual ~DistanceConstraintBond() override = default;

            /**
             * @brief Check whether a backbone bond can be formed between the two bodies.
             */
            static bool can_bond(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2);

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;
    };
}