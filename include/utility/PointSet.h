#pragma once

#include <vector>
#include <type_traits>

namespace detail {
    struct IPoint {
        static unsigned int dim() {return 0;}
    };
}

struct Point1D : detail::IPoint {
    Point1D() {}
    Point1D(double x) : x(x) {}
    Point1D(double x, double xerr) : x(x), xerr(xerr) {}

    static unsigned int dim() {return 1;}

    double x = 0, xerr = 0;
};

struct Point2D : Point1D {
    Point2D() {}
    Point2D(double x, double y) : Point1D(x), y(y) {}
    Point2D(double x, double y, double yerr) : Point1D(x), y(y), yerr(yerr) {}
    Point2D(double x, double y, double xerr, double yerr) : Point1D(x, xerr), y(y), yerr(yerr) {}

    static unsigned int dim() {return 2;}

    double y = 0, yerr = 0;
};

struct Point3D : Point2D {
    Point3D() {}
    Point3D(double x, double y, double z) : Point2D(x, y), z(z) {}

    static unsigned int dim() {return 3;}

    double z = 0, zerr = 0;
};

/**
 * @brief A representation of a set of points.
 */
template<typename T>
class PointSet {
    static_assert(std::is_base_of<detail::IPoint, T>::value);

    unsigned int dim() const noexcept {return T::dim();}

    std::vector<T> data;
};