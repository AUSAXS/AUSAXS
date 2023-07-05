#pragma once

#include <math/Vector3.h>

struct Basis3D {
    Basis3D() = default;
    Basis3D(const Vector3<double>& x, const Vector3<double>& y, const Vector3<double>& z) : x(x), y(y), z(z) {}
    Vector3<double> x, y, z;
};