// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <math/Matrix.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace ausaxs::symmetry {
    /**
     * @brief A rigid affine map  v -> rotation*v + translation.
     */
    struct AffineTransform {
        Matrix<double> rotation = Matrix<double>::identity(3);
        Vector3<double> translation{0, 0, 0};

        Vector3<double> operator()(const Vector3<double>& v) const {return rotation*v + translation;}
    };

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
     * _make_transform() to define the actual geometric mapping, allowing different symmetry
     * parameterisations without changing any of the surrounding infrastructure.
     */
    class ISymmetry {
    public:
        virtual ~ISymmetry() = default;

        /**
         * @brief Get the transform generating copy @p rep from the body's current coordinates.
         */
        [[nodiscard]] AffineTransform _get_transform(const Vector3<double>& cm, int rep = 1) const;

        /**
         * @brief Get the transform generating copy @p rep from the body's current coordinates, for a body whose orientation has
         *        changed since the symmetry was defined.
         *
         * @param body_orientation The rotation taking the owning body from the orientation it had when the symmetry was defined
         *                         to its current one. An empty optional means the body has not been reoriented.
         */
        [[nodiscard]] AffineTransform _get_transform(
            const Vector3<double>& cm, const std::optional<Matrix<double>>& body_orientation, int rep = 1
        ) const;

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
         */
        virtual std::vector<SymmetricDuplicatePair> internal_pair_schedule() const;

    protected:
        /**
         * @brief Build the map generating copy @p rep in the frame the symmetry parameters were defined in.
         *        This is the only geometry a symmetry has to supply; the public accessors above assemble the final transform from it.
         *
         * @param anchor The point to anchor the copies to, already resolved through _transform_anchor() by the caller.
         */
        [[nodiscard]] virtual AffineTransform _make_transform(const Vector3<double>& anchor, int rep) const = 0;

        /**
         * @brief The point _make_transform anchors its copies to. Defaults to the body's centre of mass; symmetries anchored
         *        elsewhere override this so the reorientation above pivots on the same point their copies do.
         *
         * Resolved once per transform and handed to _make_transform, so an override that is expensive to evaluate is not paid for twice.
         */
        [[nodiscard]] virtual Vector3<double> _transform_anchor(const Vector3<double>& cm) const;

        /**
         * @brief The orientation the copies of this symmetry are generated in, given the orientation of the body asking for them.
         *
         * Defaults to the asking body's own orientation. A symmetry shared between several bodies must instead pin every
         * participant to one common orientation, since copies placed by per-body operators would no longer form a symmetric assembly.
         */
        [[nodiscard]] virtual std::optional<Matrix<double>> _transform_orientation(const std::optional<Matrix<double>>& body_orientation) const;
    };
}
