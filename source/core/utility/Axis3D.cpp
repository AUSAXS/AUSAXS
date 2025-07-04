// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/Axis3D.h>
#include <utility/Exceptions.h>
#include <utility/Limit3D.h>
#include <math/Vector3.h>

using namespace ausaxs;

Axis3D::Axis3D() noexcept {}

Axis3D::Axis3D(const Axis3D& axis) noexcept : x(axis.x), y(axis.y), z(axis.z) {}

Axis3D::Axis3D(const Axis& x, const Axis& y, const Axis& z) noexcept : x(x), y(y), z(z) {}

Axis3D::Axis3D(const Limit3D& limits, double width) noexcept : x(limits.x, limits.x.span()/width), y(limits.y, limits.y.span()/width), z(limits.z, limits.z.span()/width) {}

Axis3D::Axis3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double width) noexcept : x(xmin, xmax, (xmax-xmin)/width), y(ymin, ymax, (ymax-ymin)/width), z(zmin, zmax, (zmax-zmin)/width) {}
Axis3D::Axis3D(const Vector3<double>& min, const Vector3<double>& max, double width) noexcept : x(min[0], max[0], (max[0]-min[0])/width), y(min[1], max[1], (max[1]-min[1])/width), z(min[2], max[2], (max[2]-min[2])/width) {}

Axis3D& Axis3D::operator=(const Axis3D& rhs) noexcept {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    return *this;
}

bool Axis3D::operator==(const Axis3D& rhs) const noexcept {
    if (x != rhs.x) {return false;}
    if (y != rhs.y) {return false;}
    if (z != rhs.z) {return false;}
    return true;
}

bool Axis3D::operator!=(const Axis3D& rhs) const noexcept {return !operator==(rhs);}

std::string Axis3D::to_string() const noexcept {
    return "Axes: \n\tx: " + x.to_string() + "\n\ty: " + y.to_string() + "\n\tz: " + z.to_string(); 
}

bool Axis3D::empty() const noexcept {return x.empty() || y.empty() || z.empty();}

void Axis3D::rebin(double width) noexcept {
    x.bins = (x.max - x.min)/width;
    y.bins = (y.max - y.min)/width;
    z.bins = (z.max - z.min)/width;
}

double Axis3D::width() const {
    if (x.bins == 0 || y.bins == 0 || z.bins == 0) {
        throw except::invalid_operation("Axis3D::width(): Cannot calculate width of an empty axis.");
    }
    auto xw = x.width();
    auto yw = y.width();
    auto zw = z.width();
    if (xw != yw || xw != zw) {
        throw except::invalid_operation("Axis3D::width(): Cannot calculate width of an axis with non-uniform bin widths.");
    }
    return xw;
}