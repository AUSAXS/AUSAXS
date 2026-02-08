// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/IDistanceConstraint.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <utility>

namespace ausaxs::rigidbody::constraints {
    class DistanceConstraintCM : public IDistanceConstraint {
        public:
            /**
             * @brief Create a center-mass constraint between a pair of bodies.
             * 
             * Complexity: O(n)
             */
            DistanceConstraintCM(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2, std::pair<int, int> isym1 = {-1, -1}, std::pair<int, int> isym2 = {-1, -1});
            virtual ~DistanceConstraintCM() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

        protected:
            double evaluate_distance() const;

        private: 
            /**
             * @brief Transforms a distance into a proper constraint for least-squares fitting. 
             * 
             * @param offset The radial offset between the new and original positions. 
             */
            static double transform(double offset);
    };
}