// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <utility/observer_ptr.h>

#include <functional>
#include <memory>

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

        struct _Relation {
            _Relation() : translation{0, 0, 0} {}
            explicit _Relation(const Vector3<double>& t) : translation(t) {}
            Vector3<double> translation;
            bool operator==(const _Relation&) const = default;
        };

        struct _Repeat {
            _Repeat() : translation{0, 0, 0}, axis{0, 0, 1}, angle(0) {}
            _Repeat(const Vector3<double>& t, const Vector3<double>& ax, double ang)
                : translation(t), axis(ax), angle(ang) {}
            _Repeat(const Vector3<double>& ax, double ang)
                : translation{0, 0, 0}, axis(ax), angle(ang) {}
            Vector3<double> translation;
            Vector3<double> axis;
            double          angle = 0;
            bool operator==(const _Repeat&) const = default;
        };

        _Relation initial_relation;
        _Repeat   repeat_relation;
        int       repetitions = 1;

        virtual std::function<Vector3<double>(Vector3<double>)> get_transform(const Vector3<double>& cm, int rep = 1) const = 0;
        virtual ISymmetry& add(observer_ptr<const ISymmetry> other) = 0;

        /**
         * @brief True if the (repetitions+1)-th copy coincides with the original body.
         */
        virtual bool is_closed() const = 0;

        virtual std::unique_ptr<ISymmetry> clone() const = 0;

        bool operator==(const ISymmetry& rhs) const = default;
    };
}
