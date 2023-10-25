#pragma once

#include <vector>
#include <iostream>

namespace detail {
    struct IPoint {
        static unsigned int dim();
    };
}

struct Point1D : ::detail::IPoint {
    Point1D();
    Point1D(double x);
    Point1D(double x, double xerr);

    static unsigned int dim();

    bool operator==(const Point1D& other) const;
    bool operator!=(const Point1D& other) const;

    double x = 0, xerr = 0;
};
std::ostream& operator<<(std::ostream& os, const Point1D& p);

struct Point2D : Point1D {
    Point2D();
    Point2D(double x, double y);
    Point2D(double x, double y, double yerr);
    Point2D(double x, double y, double xerr, double yerr);

    static unsigned int dim();

    bool operator==(const Point2D& other) const;
    bool operator!=(const Point2D& other) const;

    double y = 0, yerr = 0;
};
std::ostream& operator<<(std::ostream& os, const Point2D& p);

struct Point3D : Point2D {
    Point3D();
    Point3D(double x, double y, double z);

    static unsigned int dim();

    bool operator==(const Point3D& other) const;
    bool operator!=(const Point3D& other) const;

    double z = 0, zerr = 0;
};
std::ostream& operator<<(std::ostream& os, const Point3D& p);

/**
 * @brief A representation of a set of points.
 */
template<typename T>
class PointSet {
    static_assert(std::is_base_of<::detail::IPoint, T>::value);

    unsigned int dim() const noexcept {return T::dim();}

    std::vector<T> data;
};