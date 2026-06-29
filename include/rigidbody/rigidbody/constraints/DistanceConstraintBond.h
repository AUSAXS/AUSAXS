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

            virtual ~DistanceConstraintBond() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;
    };
}