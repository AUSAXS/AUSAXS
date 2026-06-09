// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/CyclicSymmetry.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace ausaxs::symmetry {
    /**
     * @brief A symmetry shared by several bodies: a cyclic symmetry that replicates the whole group of participating bodies as one rigid assembly.
     *
     * The rotation centre is the combined centre of mass of the participating bodies.
     *
     * The object is owned by one designated primary body; the other participating bodies must hold a ReferenceSymmetryView that forwards to it, 
     * so all of them share the same parameters.
     */
    struct ReferenceSymmetry : public ISymmetry {
        /**
         * @param bodies The participating body indices; the first is the primary that owns this symmetry.
         * @param slots  The symmetry-storage slot this symmetry occupies on each participating body (the owning slot on the primary, the view slot on the others).
         */
        ReferenceSymmetry(CyclicSymmetry base, std::vector<int> bodies, std::vector<int> slots, observer_ptr<const data::Molecule> molecule);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        std::unique_ptr<ISymmetry> clone() const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;
        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;
        std::vector<SymmetricDuplicatePair> internal_pair_schedule() const override;

        /**
         * @brief Combined centre of mass of all participating bodies (atom-count weighted).
         */
        Vector3<double> combined_cm() const;

        CyclicSymmetry base;                            //< the underlying cyclic symmetry parameters
        std::vector<int> bodies;                        //< indices of the participating bodies (primary first)
        std::vector<int> slots;                         //< symmetry slot this symmetry occupies on each body (parallel to bodies)
        observer_ptr<const data::Molecule> molecule;    //< source for the combined centre of mass
    };

    /**
     * @brief Non-owning facade so that every participating body's SymmetryStorage exposes the shared ReferenceSymmetry through the normal 
     *        ISymmetry interface. All calls are forwarded to the single owned ReferenceSymmetry.
     */
    struct ReferenceSymmetryView : public ISymmetry {
        ReferenceSymmetryView(observer_ptr<const data::Molecule> molecule, int primary_body, int symmetry_index);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        std::unique_ptr<ISymmetry> clone() const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;
        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;
        std::vector<SymmetricDuplicatePair> internal_pair_schedule() const override;

        /**
         * @brief Resolve the shared ReferenceSymmetry through the molecule. Done lazily on every call so the view survives possible body reallocations. 
         */
        observer_ptr<const ReferenceSymmetry> target() const;

        observer_ptr<const data::Molecule> molecule;    //< source for the shared symmetry
        int primary_body;                               //< body that owns the shared symmetry
        int symmetry_index;                             //< slot of the shared symmetry on the primary body
    };
}
