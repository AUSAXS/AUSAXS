#pragma once

#include <initializer_list>
#include "Vector.h"

class Vector3 : public Vector {
public:
    Vector3(const Vector3& v) : Vector(v.data) {} // copy constructor
    Vector3() : Vector(3) {} // default empty constructor
    Vector3(const std::initializer_list<double> l) : Vector(l) {} // initializer list {a, b, c, d}
    Vector3(const Vector& v) : Vector(v) {}
    Vector3(double x, double y, double z) : Vector({x, y, z}) {}
    ~Vector3() override {}

    Vector3& operator=(const Vector3& v) {
        _data = v.data;
        return *this;
    }

    Vector3& operator=(std::initializer_list<double> l) {
        _N = l.size();
        _data.assign(l);
        return *this;
    }

    double distance(const Vector3& v) const {return sqrt(pow(x-v.x, 2) + pow(y-v.y, 2) + pow(z-v.z, 2));}

    // Approximate equality, w ~ v
    bool operator==(const Vector3& v) const {return x-v.x + y-v.y + z-v.z < precision; }

    // Two Vector3s are always compatible
    void compatibility_check(const Vector3&) const {}

    double& x = _data[0];
    double& y = _data[1];
    double& z = _data[2];
};