// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/Constraint.h>
#include <data/DataFwd.h>
#include <math/Vector3.h>
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

        /**
         * @brief Capture the offset from each representative atom to the centre of mass of its body.
         *
         * A body's symmetry copies are anchored on its centre of mass, so evaluating a constraint against one needs that centre - and recomputing it on every
         * evaluation would mean a full pass over the body. The offset from the representative atom is captured once here instead: the atom and the centre both
         * move with the body, so the offset only ever rotates, and cm1()/cm2() recover the centre by turning it by however far the body has rotated since.
         *
         * Must be called once the atom indices are known and before the first distance evaluation. Assumes the body's atom composition does not change afterwards.
         */
        void cache_cm_offsets();

        /**
         * @brief Centre of mass of the constrained bodies, recovered from their representative atoms.
         */
        [[nodiscard]] Vector3<double> cm1() const;
        [[nodiscard]] Vector3<double> cm2() const; //< @copydoc cm1()

        double d_target = 0;                         // The target distance between the two bodies.
        observer_ptr<const data::Molecule> molecule; // The molecule this constraint belongs to.
        int ibody1 = -1;                             // The index of the first body.
        int ibody2 = -1;                             // The index of the second body.
        int iatom1 = -1;                             // The index of the first atom representing the CM of body1.
        int iatom2 = -1;                             // The index of the second atom representing the CM of body2.
        std::pair<int, int> isym1  = {-1, -1};       // The symmetry index of body1.
        std::pair<int, int> isym2  = {-1, -1};       // The symmetry index of body2.

        // Offset from each representative atom to its body's centre of mass, in that body's own frame. See cache_cm_offsets().
        Vector3<double> cm_offset1 = {0, 0, 0};
        Vector3<double> cm_offset2 = {0, 0, 0};
    };
}