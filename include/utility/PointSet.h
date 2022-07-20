#pragma once

#include <vector>
#include <type_traits>

namespace detail {
    class IPoint {
        public:
            static constexpr unsigned int dim() noexcept {return 0;}
    };

    class IPointError : virtual public detail::IPoint {};
}

class Point1D : virtual public detail::IPoint {
    public:
        Point1D() = default;
        Point1D(double x) : x(x) {}
        ~Point1D() = default;

        static constexpr unsigned int dim() {return 1;}

        double x = 0;
};

class Point2D : virtual public Point1D {
    public:
        Point2D() = default;
        Point2D(double x, double y) : Point1D(x), y(y) {}
        ~Point2D() = default;

        static constexpr unsigned int dim() {return 2;}

        double y = 0;
};

class Point3D : virtual public Point2D {
    public:
        Point3D() = default;
        Point3D(double x, double y, double z) : Point2D(x, y), z(z) {}
        ~Point3D() = default;

        static constexpr unsigned int dim() {return 3;}

        double z = 0;
};

class ErrorPoint1D : virtual public Point1D, virtual public detail::IPointError {
    public:
        ErrorPoint1D() = default;
        ErrorPoint1D(double x, double xerr) : Point1D(x), xerr(xerr) {}
        ~ErrorPoint1D() = default;

        static constexpr unsigned int dim() {return 1;}

        double xerr = 0;
};

class ErrorPoint2D : virtual public Point2D, virtual public ErrorPoint1D {
    public:
        ErrorPoint2D() = default;
        ErrorPoint2D(double x, double y, double xerr, double yerr) : Point2D(x, y), ErrorPoint1D(x, xerr), yerr(yerr) {}
        ~ErrorPoint2D() = default;

        static constexpr unsigned int dim() {return 2;}

        double yerr = 0;
};

class ErrorPoint3D : virtual public Point3D, virtual public ErrorPoint2D {
    public:
        ErrorPoint3D() = default;
        ErrorPoint3D(double x, double y, double z, double xerr, double yerr, double zerr) : Point3D(x, y, z), ErrorPoint2D(x, y, xerr, yerr), zerr(zerr) {}
        ~ErrorPoint3D() = default;

        static constexpr unsigned int dim() {return 3;}

        double zerr = 0;
};