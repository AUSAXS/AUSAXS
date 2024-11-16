/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <dataset/PointSet.h>

using namespace ausaxs;

unsigned int detail::IPoint::dim() {return 0;}

Point1D::Point1D() = default;
Point1D::Point1D(double x) : x(x) {}
Point1D::Point1D(double x, double xerr) : x(x), xerr(xerr) {}
unsigned int Point1D::dim() {return 1;}
bool Point1D::operator==(const Point1D& other) const {return x == other.x && xerr == other.xerr;}
bool Point1D::operator!=(const Point1D& other) const {return !(*this == other);}

Point2D::Point2D() = default;
Point2D::Point2D(double x, double y) : Point1D(x), y(y) {}
Point2D::Point2D(double x, double y, double yerr) : Point1D(x), y(y), yerr(yerr) {}
Point2D::Point2D(double x, double y, double xerr, double yerr) : Point1D(x, xerr), y(y), yerr(yerr) {}
unsigned int Point2D::dim() {return 2;}
bool Point2D::operator==(const Point2D& other) const {return Point1D::operator==(other) && y == other.y && yerr == other.yerr;}
bool Point2D::operator!=(const Point2D& other) const {return !(*this == other);}

Point3D::Point3D() = default;
Point3D::Point3D(double x, double y, double z) : Point2D(x, y), z(z) {}
unsigned int Point3D::dim() {return 3;}
bool Point3D::operator==(const Point3D& other) const {return Point2D::operator==(other) && z == other.z && zerr == other.zerr;}
bool Point3D::operator!=(const Point3D& other) const {return !(*this == other);}

std::ostream& ausaxs::operator<<(std::ostream& os, const Point1D& p) {
    os << "(" << p.x << ", " << p.xerr << ")";
    return os;
}

std::ostream& ausaxs::operator<<(std::ostream& os, const Point2D& p) {
    os << "(" << p.x << ", " << p.xerr << ", " << p.y << ", " << p.yerr << ")";
    return os;
}

std::ostream& ausaxs::operator<<(std::ostream& os, const Point3D& p) {
    os << "(" << p.x << ", " << p.xerr << ", " << p.y << ", " << p.yerr << ", " << p.z << ", " << p.zerr << ")";
    return os;
}