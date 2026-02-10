// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/Constraint.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <utility>

namespace ausaxs::rigidbody::constraints {
    struct IDistanceConstraint : Constraint {
        IDistanceConstraint() = default;
        IDistanceConstraint(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2, std::pair<int, int> isym1 = {-1, -1}, std::pair<int, int> isym2 = {-1, -1});
        IDistanceConstraint(
            observer_ptr<const data::Molecule> molecule, 
            int ibody1, int iatom1, 
            int ibody2, int iatom2,
            std::pair<int, int> isym1 = {-1, -1}, std::pair<int, int> isym2 = {-1, -1}
        );
        virtual ~IDistanceConstraint() = default;

        /**
         * @brief Get the first atom of this constraint. 
         */
        const data::AtomFF& get_atom1() const;

        /**
         * @brief Get the second atom of this constraint. 
         */
        const data::AtomFF& get_atom2() const;

        /**
         * @brief Get the first body of this constraint. 
         */
        const data::Body& get_body1() const;

        /**
         * @brief Get the second body of this constraint. 
         */
        const data::Body& get_body2() const;

        double d_target = 0;                         // The target distance between the two bodies.
        observer_ptr<const data::Molecule> molecule; // The molecule this constraint belongs to.
        int ibody1 = -1;                             // The index of the first body.
        int ibody2 = -1;                             // The index of the second body.
        int iatom1 = -1;                             // The index of the first atom representing the CM of body1.
        int iatom2 = -1;                             // The index of the second atom representing the CM of body2.
        std::pair<int, int> isym1  = {-1, -1};       // The symmetry index of body1.
        std::pair<int, int> isym2  = {-1, -1};       // The symmetry index of body2.
    };
}