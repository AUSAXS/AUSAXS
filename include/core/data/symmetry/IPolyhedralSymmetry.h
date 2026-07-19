// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/Vector3.h>
#include <math/Matrix.h>

namespace ausaxs::symmetry {
    /**
     * @brief Chiral polyhedral point-group symmetry (tetrahedral / octahedral / icosahedral).
     */
    struct IPolyhedralSymmetry : public ISymmetry {
        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;

        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;
        std::vector<SymmetricDuplicatePair> internal_pair_schedule() const override;

        Vector3<double> translation{0, 0, 0};   //< offset of the body from the group centre
        Vector3<double> rotation{0, 0, 0};      //< orientation of the group frame (Euler angles)

        /**
         * @brief The fixed rotation matrices of the point group (element 0 is the identity).
         *
         * Exposes the otherwise-protected group data so external code (e.g. the symmetry fitter)
         * can recover the group frame by matching measured rotations to these canonical elements.
         */
        [[nodiscard]] const std::vector<Matrix<double>>& group_elements() const {return group().elements;}

    protected:
        //< The fixed rotation matrices of the group (element 0 = identity) and the distance-reuse
        //< schedule derived from them. Both are invariant data, supplied once by each concrete group.
        struct GroupData {
            std::vector<Matrix<double>> elements;
            std::vector<SymmetricDuplicatePair> schedule;
        };

        virtual const GroupData& group() const = 0;

        // Build the group data from a set of generators by closing the group and bucketing its pairs.
        static GroupData build(const std::vector<Matrix<double>>& generators, int order);

        // World-space centre that every copy rotates about. The default places the body at cm + translation (offset in world coordinates). 
        virtual Vector3<double> group_centre(const Vector3<double>& cm, const Matrix<double>& F) const;
    };
}
