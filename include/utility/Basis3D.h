#pragma once

#include <math/Vector3.h>

struct Basis3D {
    Basis3D() = default;
    Basis3D(Vector3<double> x, Vector3<double> y, Vector3<double> z) : x(x), y(y), z(z) {}
    Vector3<double> x, y, z;
};