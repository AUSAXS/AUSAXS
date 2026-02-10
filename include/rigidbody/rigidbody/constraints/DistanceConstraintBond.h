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
             * @brief Create a new constraint between a pair of atoms in the two bodies.
             * 
             * Complexity: O(n)
             */
            DistanceConstraintBond(
                observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2,
                std::pair<int, int> isym1 = {-1, -1}, std::pair<int, int> isym2 = {-1, -1}
            );

            virtual ~DistanceConstraintBond() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;
    };
}