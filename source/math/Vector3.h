#pragma once

#include <initializer_list>
#include "math/Vector.h"

class Vector3 : public Vector {
public:
    using Vector::Vector;
    Vector3() : Vector(3) {}
    Vector3(Vector v) : Vector(v) {}
    Vector3(double x, double y, double z) : Vector({x, y, z}) {}

    // Assignment operator, w = v
    Vector3& operator=(const Vector3& v) {
        data.assign(v.begin(), v.end());
        return *this;
    }

    // Two Vector3s are always compatible
    void compatibility_check(const Vector3& v) const {}

    double& x = data[0];
    double& y = data[1];
    double& z = data[2];
};