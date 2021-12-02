#pragma once

#include <initializer_list>
#include "math/Vector.h"

class Vector3 : public Vector {
public:
    using Vector::Vector;
    Vector3() : Vector(3) {}
    Vector3(Vector v) : Vector(v) {}
    Vector3(double x, double y, double z) : Vector({x, y, z}) {}
    ~Vector3() override {}

    Vector3& operator=(const Vector3& v) {
        _data.assign(v.begin(), v.end());
        return *this;
    }

    Vector3& operator=(std::initializer_list<double> l) {
        _N = l.size();
        _data.assign(l);
        return *this;
    }

    // Approximate equality, w ~ v
    bool operator==(const Vector3& v) const {return x-v.x + y-v.y + z-v.z < precision; }

    // Two Vector3s are always compatible
    void compatibility_check(const Vector3& v) const {}

    double& x = _data[0];
    double& y = _data[1];
    double& z = _data[2];
};