// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/CyclicSymmetry.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace ausaxs::symmetry {
    /**
     * @brief A symmetry shared by several bodies: a cyclic symmetry that replicates the whole
     *        group of participating bodies as one rigid assembly.
     *
     * Useful when a molecule is split across multiple files (one body per file) and a single
     * symmetry should apply to all of them at once. The rotation centre is the combined
     * centre of mass of the participating bodies, so it is independent of how many bodies
     * take part — only their positions matter.
     *
     * The object is owned by one designated primary body; the other participating bodies hold
     * a ReferenceSymmetryView that forwards to it, so all of them share the same parameters.
     */
    struct ReferenceSymmetry : public ISymmetry {
        ReferenceSymmetry(CyclicSymmetry base, std::vector<int> bodies, observer_ptr<const data::Molecule> molecule);

        ISymmetry& add(observer_ptr<const ISymmetry> other) override;
        std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const override;
        std::unique_ptr<ISymmetry> clone() const override;
        unsigned int repetitions() const override;
        bool is_closed() const override;
        std::span<double> span_translation() override;
        std::span<double> span_rotation() override;
        std::vector<CopyPair> internal_pair_schedule() const override;

        /**
         * @brief Combined centre of mass of all participating bodies (atom-count weighted).
         */
        Vector3<double> combined_cm() const;

        CyclicSymmetry base;                            //< the underlying cyclic symmetry parameters
        std::vector<int> bodies;                        //< indices of the participating bodies
        observer_ptr<const data::Molecule> molecule;    //< source for the combined centre of mass
    };

    /**
     * @brief Non-owning facade so that every participating body's SymmetryStorage exposes the
     *        shared ReferenceSymmetry through the normal ISymmetry interface. All calls are
     *        forwarded to the single owned ReferenceSymmetry.
     *
     * The shared symmetry is located by (primary body index, symmetry slot) and resolved through
     * the molecule on every call rather than cached as a raw pointer. This is deliberate: the
     * rigid-body transform path replaces a body's symmetry objects (the body is restored from a
     * clone), so a cached pointer to the owning ReferenceSymmetry would dangle after the primary
     * body is transformed. The molecule and its bodies are stable for the program lifetime, so
     * resolving through them is always valid.
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
        std::vector<CopyPair> internal_pair_schedule() const override;

        /**
         * @brief Resolve the shared ReferenceSymmetry through the molecule. Done lazily on every
         *        call so the view survives body reallocation during refinement.
         */
        observer_ptr<const ReferenceSymmetry> target() const;

        observer_ptr<const data::Molecule> molecule;    //< source for the shared symmetry
        int primary_body;                               //< body that owns the shared symmetry
        int symmetry_index;                             //< slot of the shared symmetry on the primary body
    };
}
