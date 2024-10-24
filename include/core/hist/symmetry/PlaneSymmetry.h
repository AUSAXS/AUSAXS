#pragma once

#include <hist/symmetry/Symmetry.h>
#include <math/Vector3.h>

namespace ausaxs::hist::detail {class CompactCoordinates;}
namespace ausaxs::hist::symmetry {
    namespace detail {
        struct Plane {
            Plane(Vector3<double>&& x, Vector3<double>&& y) : x(std::move(x)), y(std::move(y)) {normal = x.cross(y).normalize();}
            Plane(const Vector3<double>& x, const Vector3<double>& y) : x(x), y(y) {normal = x.cross(y).normalize();}

            Vector3<double> x, y, normal;
        };
    }

    class PlaneSymmetry : public Symmetry {
        public:
            PlaneSymmetry(Vector3<double>&& x, Vector3<double>&& y) : plane(std::move(x), std::move(y)) {}
            PlaneSymmetry(const Vector3<double>& x, const Vector3<double>& y) : plane(x, y) {}

            ~PlaneSymmetry() override = default;

            /**
             * @brief Mirror a vector across the plane.
             */
            Vector3<double> mirror(const Vector3<double>& v) const;

            std::vector<double> calculate_cross_terms(hist::detail::CompactCoordinates& body1, hist::detail::CompactCoordinates& body2);

        private:
            detail::Plane plane;
    };
}