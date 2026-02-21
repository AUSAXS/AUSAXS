// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <utility/observer_ptr.h>

#include <functional>
#include <memory>
#include <span>

namespace ausaxs::symmetry {
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
        bool operator==(const ISymmetry& rhs) const = default;

        virtual std::span<double> span_translation() = 0;
        virtual std::span<double> span_rotation() = 0;
    };
}
