#include "math/Vector3.h"
#include "math/Matrix.h"

Vector3::Vector3(Vector3&& v) noexcept : Vector(std::move(v)) {} // move constructor
Vector3::Vector3(const Vector3& v) : Vector(v.data) {} // copy constructor
Vector3::Vector3() : Vector(3) {} // default empty constructor
Vector3::Vector3(const std::initializer_list<double> l) : Vector(l) {} // initializer list {a, b, c, d}
Vector3::Vector3(const Vector& v) : Vector(v) {}
Vector3::Vector3(double x, double y, double z) : Vector({x, y, z}) {}
Vector3::~Vector3() = default;

Vector3& Vector3::operator=(const Vector3& v) {
    data = v.data;
    return *this;
}

Vector3& Vector3::operator=(std::initializer_list<double> l) {
    N = l.size();
    data.assign(l);
    return *this;
}

double Vector3::distance(const Vector3& v) const {return sqrt(pow(x()-v.x(), 2) + pow(y()-v.y(), 2) + pow(z()-v.z(), 2));}

Vector3 Vector3::cross(const Vector3& v) const {return {y()*v.z() - v.y()*z(), z()*v.x() - v.z()*x(), x()*v.y() - v.x()*y()};}

bool Vector3::operator==(const Vector3& v) const {return abs(x()-v.x()) + abs(y()-v.y()) + abs(z()-v.z()) < precision; }

void Vector3::compatibility_check(const Vector3&) const {}

Vector3 Vector3::normalize() {
    operator=(operator/(norm()));
    return *this;
}

Vector3 Vector3::normalize() const {
    return operator/(norm());
}

Vector3 Vector3::normalize_copy() const {
    return operator/(norm());
}

std::tuple<Vector3, Vector3, Vector3> Vector3::generate_basis(const Vector3& v) {
    Vector3 n = v.normalize_copy();

    // Handle the singularity
    if (n.z() < -0.9999999) { 
        Vector3 b1(0, -1, 0);
        Vector3 b2(-1, 0, 0);
        return std::make_tuple(n, b1, b2);
    }
    const float a = 1/(1 + n.z());
    const float b = -n.x()*n.y()*a;
    Vector3 b1(1-n.x()*n.x()*a, b, -n.x());
    Vector3 b2(b, 1-n.y() * n.y()*a, -n.y());
    return std::make_tuple(n, b1, b2);
}

std::tuple<Vector3, Vector3, Vector3> Vector3::generate_basis() {
    return Vector3::generate_basis(*this);
}

void Vector3::rotate(const Matrix<double>& matrix) {
    *this = matrix*(*this);
}

void Vector3::rotate(const Vector3& axis, double angle) {
    Matrix R = Matrix<double>::rotation_matrix(axis, angle);
    rotate(R);
}