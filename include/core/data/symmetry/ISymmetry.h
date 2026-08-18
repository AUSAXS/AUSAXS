// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <utility/observer_ptr.h>

#include <functional>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace ausaxs::symmetry {
    /**
     * @brief One representative inter-copy distance correlation job.
     *
     * repA / repB are repetition indices: 0 denotes the original body, 1..repetitions()
     * denote the generated copies. scale is the multiplicity: how many identical copy-pairs
     * this single representative stands for. Every copy-pair not listed in a schedule is
     * geometrically identical to one that is.
     */
    struct SymmetricDuplicatePair {
        int repA;
        int repB;
        int scale;
    };

    /**
     * @brief Abstract base class for symmetry operations used in rigid-body optimisation.
     *
     * Holds the shared parameter storage (initial_relation, repeat_relation, repetitions)
     * that both the optimiser and histogram backend read and write.  Subclasses implement
     * get_transform() to define the actual geometric mapping, allowing different symmetry
     * parameterisations without changing any of the surrounding infrastructure.
     */
    class ISymmetry {
    public:
        virtual ~ISymmetry() = default;

        virtual std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const = 0;
        virtual ISymmetry& add(observer_ptr<const ISymmetry> other) = 0;

        /**
         * @brief True if the (repetitions+1)-th copy coincides with the original body.
         */
        virtual bool is_closed() const = 0;
        virtual unsigned int repetitions() const = 0;
        virtual std::unique_ptr<ISymmetry> clone() const = 0;
        bool operator==(const ISymmetry& rhs) const = delete;

        virtual std::span<double> span_translation() = 0;
        virtual std::span<double> span_rotation() = 0;

        /**
         * @brief Express a rotation of the given angle as a delta in this symmetry's own rotation parameters.
         *
         * Lets a sampler work in radians without knowing how each symmetry stores its rotation. The returned value is
         * meant to be written into span_rotation() and added to the live parameters; it is not applied here.
         * The default treats the rotation parameters as a rotation vector, which is what every symmetry
         * parameterised by Euler angles needs. CyclicSymmetry, whose parameter is an axis direction, overrides it.
         *
         * @param angle The rotation angle, in radians.
         * @param direction A unit vector giving the rotation direction. Supplied by the caller so that symmetries need no RNG of their own.
         */
        virtual Vector3<double> rotation_from_angle(double angle, const Vector3<double>& direction) const;

        /**
         * @brief Short predefined-name string identifying this symmetry's type (e.g. "c4", "p2", "d3"),
         *        matching the names accepted by symmetry::get(std::string_view) / symmetry::create().
         *        Composite symmetries return their nested "inner-outer" name (e.g. "p2-c3").
         */
        virtual std::string type_name() const = 0;

        /**
         * @brief Distinct inter-copy distance pairs within {original, copy_1, ..., copy_N}.
         *
         * The histogram backend evaluates one cross-correlation per returned CopyPair and
         * weights it by CopyPair::scale; every other copy-pair is identical to a listed one.
         * The default implementation reproduces the cyclic-chain reuse: copy k sits one
         * fixed generator step from copy k-1, so the distance depends only on the index
         * separation. Subclasses with a different group structure override this.
         */
        virtual std::vector<SymmetricDuplicatePair> internal_pair_schedule() const;
    };
}
